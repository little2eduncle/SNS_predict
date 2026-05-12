function [C_total, C_windows, C_intensity, C_XX_intensity] = ...
    simulate_SNS_TF_QKD_2(N, thetaA, thetaB, thetaC, varargin)
% simulate_SNS_TF_QKD (fast GPU version, no full D0/D1 output)
%
% 主要优化：
%   1) 数值计算放到 GPU
%   2) 不返回 D0/D1 全长数组
%   3) 不分别返回 A_isZ/B_isZ/send/type/idx 等多个大数组
%   4) 每个 chunk 只 gather 一个小的 uint8 打包矩阵
%
% 这版特别修正：
%   XX 窗口统计包含 [decoys..., vacuum] x [decoys..., vacuum]
%   即同时统计：
%       decoy-decoy, decoy-vacuum, vacuum-decoy, vacuum-vacuum
%
% category 编码:
%   1 = Z-signal
%   2 = Z-vacuum
%   3..(2+n_x) = X-window intensity index over x_set = [decoys..., vacuum]
%
% outcome 编码:
%   1 -> 00
%   2 -> 01 (only D1)
%   3 -> 10 (only D0)
%   4 -> 11
%
% 输入:
%   thetaA = {signal_A, decoy_A, alpha_A, dAC}
%   thetaB = {signal_B, decoy_B, alpha_B, dBC}
%   thetaC = {pa0, pa1, pc0, pc1, pd0, pd1}
%
% 可选参数:
%   'buffer'      : chunk size, default 1e4
%   'z_window'    : probability of choosing Z window
%   'epsilon'     : probability of sending signal in Z window
%   'p_decoy'     : probability vector over [decoy_1, decoy_2, ..., vacuum]
%   'lambda'      : phase-slice threshold
%   'fluct_sigma' : Gaussian std for nonzero intensities
%   'vacuum'      : vacuum intensity (default 0)
%   'use_gpu'     : true/false, default true

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'buffer', 1e4);
addParameter(p, 'z_window', 0.5);
addParameter(p, 'epsilon', 0.99);
addParameter(p, 'p_decoy', []);
addParameter(p, 'lambda', 0.2);
addParameter(p, 'fluct_sigma', 0.0);
addParameter(p, 'vacuum', 0);
addParameter(p, 'use_gpu', true);
parse(p, varargin{:});
opt = p.Results;

buffer      = opt.buffer;
z_window    = opt.z_window;
epsilon     = opt.epsilon;
p_decoy     = opt.p_decoy;
lambda      = opt.lambda;
fluct_sigma = opt.fluct_sigma;
vacuum      = opt.vacuum;
use_gpu     = opt.use_gpu;

if use_gpu
    assert(gpuDeviceCount > 0, '未检测到可用 GPU。请设置 ''use_gpu'', false。');
    gpuDevice;
end

% -------------------------------------------------------------------------
% Unpack thetaA / thetaB / thetaC
% -------------------------------------------------------------------------
signal_A = thetaA{1};
decoy_A  = thetaA{2};
alpha_A  = thetaA{3};
dAC      = thetaA{4};

signal_B = thetaB{1};
decoy_B  = thetaB{2};
alpha_B  = thetaB{3};
dBC      = thetaB{4};

pa0 = thetaC{1};
pa1 = thetaC{2};
pc0 = thetaC{3};
pc1 = thetaC{4};
pd0 = thetaC{5};
pd1 = thetaC{6};

decoy_A = decoy_A(:).';
decoy_B = decoy_B(:).';

x_set_A = [decoy_A, vacuum];
x_set_B = [decoy_B, vacuum];

if isempty(p_decoy)
    p_decoy = ones(1, numel(x_set_A)) / numel(x_set_A);
end

% -------------------------------------------------------------------------
% Checks
% -------------------------------------------------------------------------
assert(isscalar(signal_A), 'signal_A must be a scalar.');
assert(isscalar(signal_B), 'signal_B must be a scalar.');
assert(abs(sum(p_decoy) - 1) < 1e-12, 'The elements of p_decoy must sum to 1.');
assert(numel(x_set_A) == numel(p_decoy), ...
    'Length mismatch: p_decoy must match [decoy_A, vacuum].');
assert(numel(x_set_B) == numel(p_decoy), ...
    'Length mismatch: p_decoy must match [decoy_B, vacuum].');
assert(numel(x_set_A) == numel(x_set_B), ...
    'Current version assumes Alice and Bob use the same number of X-window intensities.');

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
n_signalA = 2;
n_signalB = 2;
n_xA      = numel(x_set_A);
n_xB      = numel(x_set_B);

C_total = zeros(1, 4);

C_ZZ = zeros(1, 4);
C_ZX = zeros(1, 4);
C_XZ = zeros(1, 4);
C_XX = zeros(1, 4);

C_intensity.signal_signal = zeros(n_signalA, n_signalB, 4);
C_intensity.signal_decoy  = zeros(n_signalA, n_xB, 4);
C_intensity.decoy_signal  = zeros(n_xA, n_signalB, 4);
C_intensity.decoy_decoy   = zeros(n_xA, n_xB, 4);

% 这里改成完整 X×X 组合: [decoys..., vacuum] × [decoys..., vacuum]
C_XX_intensity.mu_nominal_A = x_set_A(:);
C_XX_intensity.mu_nominal_B = x_set_B(:).';
C_XX_intensity.counts       = zeros(n_xA, n_xB, 4);
C_XX_intensity.D0_only      = zeros(n_xA, n_xB);
C_XX_intensity.D1_only      = zeros(n_xA, n_xB);
C_XX_intensity.single_click = zeros(n_xA, n_xB);
C_XX_intensity.double_click = zeros(n_xA, n_xB);
C_XX_intensity.total_valid  = zeros(n_xA, n_xB);

num_chunks = ceil(N / buffer);

etaA = 10^(-alpha_A * dAC / 10);
etaB = 10^(-alpha_B * dBC / 10);

if use_gpu
    etaA = gpuArray(etaA);
    etaB = gpuArray(etaB);

    pa0 = gpuArray(pa0); pa1 = gpuArray(pa1);
    pc0 = gpuArray(pc0); pc1 = gpuArray(pc1);
    pd0 = gpuArray(pd0); pd1 = gpuArray(pd1);
    lambda = gpuArray(lambda);

    x_set_A_dev   = gpuArray(x_set_A);
    x_set_B_dev   = gpuArray(x_set_B);
    cdf_decoy_dev = gpuArray(cumsum(p_decoy));
else
    x_set_A_dev   = x_set_A;
    x_set_B_dev   = x_set_B;
    cdf_decoy_dev = cumsum(p_decoy);
end

if exist('start_up', 'file') == 2
    start_up;
end

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------

start_up;

for chunk = 1:num_chunks
    first = (chunk - 1) * buffer + 1;
    last  = min(chunk * buffer, N);
    chunk_size = last - first + 1;

    if use_gpu
        A_isZ = rand(1, chunk_size, 'gpuArray') < z_window;
        B_isZ = rand(1, chunk_size, 'gpuArray') < z_window;
    else
        A_isZ = rand(1, chunk_size) < z_window;
        B_isZ = rand(1, chunk_size) < z_window;
    end

    % 只返回 mu 和 category code
    [muA_nominal, catA] = prepare_party_states_cat( ...
        A_isZ, signal_A, vacuum, x_set_A_dev, epsilon, cdf_decoy_dev, use_gpu);

    [muB_nominal, catB] = prepare_party_states_cat( ...
        B_isZ, signal_B, vacuum, x_set_B_dev, epsilon, cdf_decoy_dev, use_gpu);

    % 直接生成相位差
    if use_gpu
        delta_theta = 2 * pi * rand(1, chunk_size, 'gpuArray');
    else
        delta_theta = 2 * pi * rand(1, chunk_size);
    end

    muA = apply_intensity_fluctuation_local_gpu(muA_nominal, fluct_sigma);
    muB = apply_intensity_fluctuation_local_gpu(muB_nominal, fluct_sigma);

    muA_eff = muA .* etaA;
    muB_eff = muB .* etaB;

    cdt = cos(delta_theta);
    inter = 2 .* sqrt(muA_eff .* muB_eff) .* cdt;

    I0 = 0.5 .* (muA_eff + muB_eff + inter);
    I1 = 0.5 .* (muA_eff + muB_eff - inter);

    I0 = max(I0, 0);
    I1 = max(I1, 0);

    p0_final = 1 - (1 - pd0) .* (1 - pa0) .* exp(-pc0 .* I0);
    p1_final = 1 - (1 - pd1) .* (1 - pa1) .* exp(-pc1 .* I1);

    p0_final = min(max(p0_final, 0), 1);
    p1_final = min(max(p1_final, 0), 1);

    if use_gpu
        finalD0_chunk = rand(1, chunk_size, 'gpuArray') < p0_final;
        finalD1_chunk = rand(1, chunk_size, 'gpuArray') < p1_final;
    else
        finalD0_chunk = rand(1, chunk_size) < p0_final;
        finalD1_chunk = rand(1, chunk_size) < p1_final;
    end

    % 1 -> 00, 2 -> 01, 3 -> 10, 4 -> 11
    outcome = 1 + double(finalD1_chunk) + 2 * double(finalD0_chunk);

    % XX 窗口有效事件：所有 X×X 组合，只要通过 phase slice
    phase_ok = abs(cdt) >= (1 - lambda);
    valid_xx = (catA >= 3) & (catB >= 3) & phase_ok;

    % 只 gather 小打包矩阵
    packed = gather(uint8([catA; catB; outcome; valid_xx]));

    A = packed(1, :);        % category A
    B = packed(2, :);        % category B
    O = packed(3, :);        % outcome
    V = packed(4, :) > 0;    % valid XX

    % ---------------------------------------------------------------------
    % Total / window counts
    % ---------------------------------------------------------------------
    C_total = C_total + count_outcomes4(O);

    % Z-window:
    %   1 = Z-signal
    %   2 = Z-vacuum
    % 在 ZZ 中保留 xor(sendA, sendB) => (1,2) or (2,1)
    maskZZ = (A <= 2) & (B <= 2) & (A ~= B);
    maskZX = (A <= 2) & (B >= 3);
    maskXZ = (A >= 3) & (B <= 2);

    C_ZZ = C_ZZ + count_outcomes4(O(maskZZ));
    C_ZX = C_ZX + count_outcomes4(O(maskZX));
    C_XZ = C_XZ + count_outcomes4(O(maskXZ));
    C_XX = C_XX + count_outcomes4(O(V));

    % ---------------------------------------------------------------------
    % Counts by intensity combination
    % ---------------------------------------------------------------------
    maskSS = (A <= 2) & (B <= 2);
    if any(maskSS)
        ia = double(A(maskSS));      % 1..2
        ib = double(B(maskSS));      % 1..2
        oo = double(O(maskSS));
        C_intensity.signal_signal = C_intensity.signal_signal + ...
            accum3_count(ia, ib, oo, [n_signalA, n_signalB, 4]);
    end

    maskSD = (A <= 2) & (B >= 3);
    if any(maskSD)
        ia = double(A(maskSD));          % 1..2
        ib = double(B(maskSD) - 2);      % 1..n_xB
        oo = double(O(maskSD));
        C_intensity.signal_decoy = C_intensity.signal_decoy + ...
            accum3_count(ia, ib, oo, [n_signalA, n_xB, 4]);
    end

    maskDS = (A >= 3) & (B <= 2);
    if any(maskDS)
        ia = double(A(maskDS) - 2);      % 1..n_xA
        ib = double(B(maskDS));          % 1..2
        oo = double(O(maskDS));
        C_intensity.decoy_signal = C_intensity.decoy_signal + ...
            accum3_count(ia, ib, oo, [n_xA, n_signalB, 4]);
    end

    maskDD = (A >= 3) & (B >= 3);
    if any(maskDD)
        ia = double(A(maskDD) - 2);      % 1..n_xA
        ib = double(B(maskDD) - 2);      % 1..n_xB
        oo = double(O(maskDD));
        C_intensity.decoy_decoy = C_intensity.decoy_decoy + ...
            accum3_count(ia, ib, oo, [n_xA, n_xB, 4]);
    end

    % ---------------------------------------------------------------------
    % Valid XX counts for all X×X combinations, including vacuum combos
    % ---------------------------------------------------------------------
    if any(V)
        ia = double(A(V) - 2);           % 1..n_xA
        ib = double(B(V) - 2);           % 1..n_xB
        oo = double(O(V));

        C_XX_intensity.counts = C_XX_intensity.counts + ...
            accum3_count(ia, ib, oo, [n_xA, n_xB, 4]);

        C_XX_intensity.total_valid = C_XX_intensity.total_valid + ...
            accum2_count(ia, ib, [n_xA, n_xB]);

        maskD1 = (oo == 2);
        if any(maskD1)
            C_XX_intensity.D1_only = C_XX_intensity.D1_only + ...
                accum2_count(ia(maskD1), ib(maskD1), [n_xA, n_xB]);
        end

        maskD0 = (oo == 3);
        if any(maskD0)
            C_XX_intensity.D0_only = C_XX_intensity.D0_only + ...
                accum2_count(ia(maskD0), ib(maskD0), [n_xA, n_xB]);
        end

        maskSingle = (oo == 2) | (oo == 3);
        if any(maskSingle)
            C_XX_intensity.single_click = C_XX_intensity.single_click + ...
                accum2_count(ia(maskSingle), ib(maskSingle), [n_xA, n_xB]);
        end

        maskDouble = (oo == 4);
        if any(maskDouble)
            C_XX_intensity.double_click = C_XX_intensity.double_click + ...
                accum2_count(ia(maskDouble), ib(maskDouble), [n_xA, n_xB]);
        end
    end
end

C_windows = struct('ZZ', C_ZZ, ...
                   'ZX', C_ZX, ...
                   'XZ', C_XZ, ...
                   'XX', C_XX);

end

% =========================================================================
% State preparation -> category code
% cat:
%   1 = Z-signal
%   2 = Z-vacuum
%   3..(2+n_x) = X-window x_set index
% =========================================================================
function [mu, cat] = prepare_party_states_cat( ...
    isZ, signal_mu, vacuum_mu, x_set, epsilon, cdf_decoy, use_gpu)

n = numel(isZ);

if use_gpu
    mu  = zeros(1, n, 'gpuArray');
    cat = zeros(1, n, 'gpuArray');
    rZ  = rand(1, n, 'gpuArray');
    rX  = rand(1, n, 'gpuArray');
else
    mu  = zeros(1, n);
    cat = zeros(1, n);
    rZ  = rand(1, n);
    rX  = rand(1, n);
end

% Z window
zSend = isZ & (rZ < epsilon);
zVac  = isZ & ~zSend;

mu(zSend)  = signal_mu;
cat(zSend) = 1;

mu(zVac)   = vacuum_mu;
cat(zVac)  = 2;

% X window
isX = ~isZ;

lower = 0;
for k = 1:(numel(x_set) - 1)
    mask = isX & (rX > lower) & (rX <= cdf_decoy(k));
    mu(mask)  = x_set(k);
    cat(mask) = k + 2;
    lower = cdf_decoy(k);
end

% 最后一档直接吃掉剩余概率
maskLast = isX & (cat == 0);
mu(maskLast)  = x_set(end);
cat(maskLast) = numel(x_set) + 2;

end

% =========================================================================
% Gaussian fluctuation
% =========================================================================
function mu_out = apply_intensity_fluctuation_local_gpu(mu_in, fluct_sigma)

if fluct_sigma <= 0
    mu_out = mu_in;
    return;
end

mu_out = mu_in;
mask = mu_in > 0;

if isa(mu_in, 'gpuArray')
    noise = fluct_sigma .* randn(size(mu_in), 'like', mu_in);
else
    noise = fluct_sigma .* randn(size(mu_in));
end

mu_out(mask) = mu_in(mask) + noise(mask);
mu_out(mask) = max(mu_out(mask), 0);

end

% =========================================================================
% Count helpers
% =========================================================================
function c = count_outcomes4(outcome_vec)
if isempty(outcome_vec)
    c = zeros(1, 4);
else
    c = accumarray(double(outcome_vec(:)), 1, [4, 1]).';
end
end

function A = accum2_count(i1, i2, out_sz)
if isempty(i1)
    A = zeros(out_sz);
else
    A = accumarray([i1(:), i2(:)], 1, out_sz);
end
end

function A = accum3_count(i1, i2, i3, out_sz)
if isempty(i1)
    A = zeros(out_sz);
else
    A = accumarray([i1(:), i2(:), i3(:)], 1, out_sz);
end
end