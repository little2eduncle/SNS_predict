function [C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
    simulate_SNS_TF_QKD_5(N, thetaA, thetaB, thetaC, varargin)
% simulate_SNS_TF_QKD
%
% 功能:
%   1) 记录 ZZ 窗口下“所有”光强组合的三类事件 {01, 10, rest}
%   2) 记录 XX 窗口下“通过相位筛选”的所有光强组合的三类事件 {01, 10, rest}
%   3) 额外记录 XX 窗口下“筛选前”的各光强组合总出现次数 total_raw
%
% 其中:
%   01   = 只有 D1 响应
%   10   = 只有 D0 响应
%   rest = {00,11}，即“都不响应 + 都响应”合并计数
%
% -------------------------------------------------------------------------
% category 编码:
%   1 = Z-signal
%   2 = Z-vacuum
%   3..(2+n_x) = X-window intensity index over x_set = [decoys..., vacuum]
%
% outcome 编码:
%   1 -> 00   (都不响应)
%   2 -> 01   (只有 D1 响应)
%   3 -> 10   (只有 D0 响应)
%   4 -> 11   (D0 和 D1 都响应)
%
% counts3 编码:
%   counts3(:,:,1) -> 01   (only D1)
%   counts3(:,:,2) -> 10   (only D0)
%   counts3(:,:,3) -> rest = {00,11}
%
% -------------------------------------------------------------------------
% 输入:
%   N      : 总脉冲数
%
%   thetaA = {signal_A, decoy_A, alpha_A, dAC}
%       signal_A : Alice 在 Z 窗口下的 signal 光强（标量）
%       decoy_A  : Alice 在 X 窗口下的 decoy 光强数组
%       alpha_A  : Alice-Charlie 光纤损耗系数 (dB/km)
%       dAC      : Alice 到 Charlie 的距离 (km)
%
%   thetaB = {signal_B, decoy_B, alpha_B, dBC}
%       signal_B : Bob 在 Z 窗口下的 signal 光强（标量）
%       decoy_B  : Bob 在 X 窗口下的 decoy 光强数组
%       alpha_B  : Bob-Charlie 光纤损耗系数 (dB/km)
%       dBC      : Bob 到 Charlie 的距离 (km)
%
%   thetaC = {pa0, pa1, pc0, pc1, pd0, pd1}
%       pa0, pa1 : D0 / D1 的后脉冲概率
%       pc0, pc1 : D0 / D1 的有效探测系数
%       pd0, pd1 : D0 / D1 的暗计数概率
%
% -------------------------------------------------------------------------
% 可选参数:
%   'buffer'      : 每个 chunk 的大小，默认 1e4
%   'z_window'    : 选择 Z 窗口的概率
%   'epsilon'     : 在 Z 窗口中发送 signal 的概率
%   'p_decoy'     : X 窗口中 [decoy_1, decoy_2, ..., vacuum] 的选择概率
%   'lambda'      : XX 窗口相位筛选阈值，条件为 abs(cos(delta_theta)) >= 1-lambda
%   'fluct_sigma' : 非零光强的加性高斯波动标准差；
%                   具体实现为 mu_out = max(mu_in + N(0, fluct_sigma^2), 0)
%                   真空态 (mu_in = 0) 不加波动
%   'vacuum'      : 真空态光强，默认 0
%   'use_gpu'     : 是否使用 GPU，默认 true
%
% -------------------------------------------------------------------------
% 输出:
%
%   C_total : 1x4 向量
%       所有脉冲的四类 outcome 总计数:
%           [N_00, N_01, N_10, N_11]
%
%   C_windows : struct
%       窗口级别的四类 outcome 计数:
%           C_windows.ZZ = [N_00, N_01, N_10, N_11]
%               -> 所有 ZZ 组合的总计数
%           C_windows.XX = [N_00, N_01, N_10, N_11]
%               -> 所有“通过相位筛选”的 XX 组合总计数
%
%   C_ZZ_intensity : struct
%       ZZ 窗口下所有 2x2 光强组合的三类事件统计
%
%   C_XX_intensity : struct
%       XX 窗口下所有 X 光强组合的三类事件统计
%
%       字段:
%           .mu_nominal_A   : [n_xA x 1]，Alice 的 X 窗口名义光强
%           .mu_nominal_B   : [1 x n_xB]，Bob 的 X 窗口名义光强
%           .counts3        : [n_xA x n_xB x 3]
%           .D1_only        : [n_xA x n_xB]
%           .D0_only        : [n_xA x n_xB]
%           .rest           : [n_xA x n_xB]
%           .total_valid    : [n_xA x n_xB]
%                             表示该 XX 光强组合中，通过相位筛选后的总次数
%           .total_raw      : [n_xA x n_xB]
%                             表示该 XX 光强组合中，筛选前（所有 XX 事件）的总次数

% -------------------------------------------------------------------------
% 解析输入参数
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'buffer', 1e4);
addParameter(p, 'z_window', 0.8);
addParameter(p, 'epsilon', 0.05);
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

if exist('start_up', 'file') == 2
    start_up;
end

if use_gpu
    assert(gpuDeviceCount > 0, '未检测到可用 GPU。请设置 ''use_gpu'', false。');
    gpuDevice;
end

% -------------------------------------------------------------------------
% 读取参数
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
% 合法性检查
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
% 初始化统计结构
% -------------------------------------------------------------------------
n_signalA = 2;
n_signalB = 2;
n_xA      = numel(x_set_A);
n_xB      = numel(x_set_B);

C_total = zeros(1, 4);
C_ZZ = zeros(1, 4);
C_XX = zeros(1, 4);

% ZZ
C_ZZ_intensity.mu_nominal_A = [signal_A; vacuum];
C_ZZ_intensity.mu_nominal_B = [signal_B, vacuum];
C_ZZ_intensity.counts3     = zeros(n_signalA, n_signalB, 3);
C_ZZ_intensity.D1_only     = zeros(n_signalA, n_signalB);
C_ZZ_intensity.D0_only     = zeros(n_signalA, n_signalB);
C_ZZ_intensity.rest        = zeros(n_signalA, n_signalB);
C_ZZ_intensity.total_valid = zeros(n_signalA, n_signalB);

% XX
C_XX_intensity.mu_nominal_A = x_set_A(:);
C_XX_intensity.mu_nominal_B = x_set_B(:).';
C_XX_intensity.counts3     = zeros(n_xA, n_xB, 3);
C_XX_intensity.D1_only     = zeros(n_xA, n_xB);
C_XX_intensity.D0_only     = zeros(n_xA, n_xB);
C_XX_intensity.rest        = zeros(n_xA, n_xB);
C_XX_intensity.total_valid = zeros(n_xA, n_xB);
C_XX_intensity.total_raw   = zeros(n_xA, n_xB);

num_chunks = ceil(N / buffer);

% -------------------------------------------------------------------------
% 固定信道透过率
% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
% 主循环：分块仿真
% -------------------------------------------------------------------------
count = 0;
for chunk = 1:num_chunks

    first = (chunk - 1) * buffer + 1;
    last  = min(chunk * buffer, N);
    chunk_size = last - first + 1;

    % 1) 窗口选择
    if use_gpu
        A_isZ = rand(1, chunk_size, 'gpuArray') < z_window;
        B_isZ = rand(1, chunk_size, 'gpuArray') < z_window;
    else
        A_isZ = rand(1, chunk_size) < z_window;
        B_isZ = rand(1, chunk_size) < z_window;
    end

    % 2) 状态选择
    [muA_nominal, catA] = prepare_party_states_cat( ...
        A_isZ, signal_A, vacuum, x_set_A_dev, epsilon, cdf_decoy_dev, use_gpu);

    [muB_nominal, catB] = prepare_party_states_cat( ...
        B_isZ, signal_B, vacuum, x_set_B_dev, epsilon, cdf_decoy_dev, use_gpu);

    % 3) 相位差
    if use_gpu
        delta_theta = 2 * pi * rand(1, chunk_size, 'gpuArray');
    else
        delta_theta = 2 * pi * rand(1, chunk_size);
    end

    % 4) 光强波动：加性高斯 + 裁到非负；真空态不波动
    muA = apply_intensity_fluctuation_local_gpu(muA_nominal, fluct_sigma);
    muB = apply_intensity_fluctuation_local_gpu(muB_nominal, fluct_sigma);

    % 5) 经过信道衰减后的有效光强（已含探测系数）
    muA_eff_I0 = muA .* pc0 .* etaA;
    muA_eff_I1 = muA .* pc1 .* etaA;

    muB_eff_I0 = muB .* pc0 .* etaB;
    muB_eff_I1 = muB .* pc1 .* etaB;

    % 6) 50:50 BS 干涉
    cdt = cos(delta_theta);
    inter0 = 2 .* sqrt(max(muA_eff_I0 .* muB_eff_I0, 0)) .* cdt;
    inter1 = 2 .* sqrt(max(muA_eff_I1 .* muB_eff_I1, 0)) .* cdt;

    I0 = 0.5 .* (muA_eff_I0 + muB_eff_I0 + inter0);
    I1 = 0.5 .* (muA_eff_I1 + muB_eff_I1 - inter1);

    I0 = max(I0, 0);
    I1 = max(I1, 0);

    count = count + sum(muA_nominal, 2) + sum(muB_nominal, 2);

    % 7) 点击概率
    p0_final = 1 - (1 - pd0) .* (1 - pa0) .* exp(-I0);
    p1_final = 1 - (1 - pd1) .* (1 - pa1) .* exp(-I1);

    p0_final = min(max(p0_final, 0), 1);
    p1_final = min(max(p1_final, 0), 1);

    % 8) 伯努利抽样
    if use_gpu
        finalD0_chunk = rand(1, chunk_size, 'gpuArray') < p0_final;
        finalD1_chunk = rand(1, chunk_size, 'gpuArray') < p1_final;
    else
        finalD0_chunk = rand(1, chunk_size) < p0_final;
        finalD1_chunk = rand(1, chunk_size) < p1_final;
    end

    outcome = 1 + double(finalD1_chunk) + 2 * double(finalD0_chunk);

    % 9) ZZ / XX 统计口径
    phase_ok = abs(cdt) >= (1 - lambda);
    raw_xx   = (catA >= 3) & (catB >= 3);
    valid_xx = raw_xx & phase_ok;
    valid_zz = (catA <= 2) & (catB <= 2);

    % 10) gather 小矩阵
    packed = gather(uint8([catA; catB; outcome; valid_zz; valid_xx; raw_xx]));

    A    = packed(1, :);
    Bc   = packed(2, :);
    O    = packed(3, :);
    VZZ  = packed(4, :) > 0;
    VXX  = packed(5, :) > 0;
    VRAW = packed(6, :) > 0;

    % 11) 总体 / 窗口统计
    C_total = C_total + count_outcomes4(O);
    C_ZZ    = C_ZZ + count_outcomes4(O(VZZ));
    C_XX    = C_XX + count_outcomes4(O(VXX));

    % 11.5) XX 筛选前原始总数
    if any(VRAW)
        ia_raw = double(A(VRAW)  - 2);
        ib_raw = double(Bc(VRAW) - 2);
        C_XX_intensity.total_raw = C_XX_intensity.total_raw + ...
            accum2_count(ia_raw, ib_raw, [n_xA, n_xB]);
    end

    % 12) ZZ 统计
    if any(VZZ)
        ia = double(A(VZZ));
        ib = double(Bc(VZZ));
        oo = double(O(VZZ));

        C_ZZ_intensity.total_valid = C_ZZ_intensity.total_valid + ...
            accum2_count(ia, ib, [n_signalA, n_signalB]);

        maskD1 = (oo == 2);
        if any(maskD1)
            tmp = accum2_count(ia(maskD1), ib(maskD1), [n_signalA, n_signalB]);
            C_ZZ_intensity.D1_only        = C_ZZ_intensity.D1_only + tmp;
            C_ZZ_intensity.counts3(:,:,1) = C_ZZ_intensity.counts3(:,:,1) + tmp;
        end

        maskD0 = (oo == 3);
        if any(maskD0)
            tmp = accum2_count(ia(maskD0), ib(maskD0), [n_signalA, n_signalB]);
            C_ZZ_intensity.D0_only        = C_ZZ_intensity.D0_only + tmp;
            C_ZZ_intensity.counts3(:,:,2) = C_ZZ_intensity.counts3(:,:,2) + tmp;
        end

        maskRest = (oo == 1) | (oo == 4);
        if any(maskRest)
            tmp = accum2_count(ia(maskRest), ib(maskRest), [n_signalA, n_signalB]);
            C_ZZ_intensity.rest           = C_ZZ_intensity.rest + tmp;
            C_ZZ_intensity.counts3(:,:,3) = C_ZZ_intensity.counts3(:,:,3) + tmp;
        end
    end

    % 13) XX 统计（仅通过筛选）
    if any(VXX)
        ia = double(A(VXX)  - 2);
        ib = double(Bc(VXX) - 2);
        oo = double(O(VXX));

        C_XX_intensity.total_valid = C_XX_intensity.total_valid + ...
            accum2_count(ia, ib, [n_xA, n_xB]);

        maskD1 = (oo == 2);
        if any(maskD1)
            tmp = accum2_count(ia(maskD1), ib(maskD1), [n_xA, n_xB]);
            C_XX_intensity.D1_only        = C_XX_intensity.D1_only + tmp;
            C_XX_intensity.counts3(:,:,1) = C_XX_intensity.counts3(:,:,1) + tmp;
        end

        maskD0 = (oo == 3);
        if any(maskD0)
            tmp = accum2_count(ia(maskD0), ib(maskD0), [n_xA, n_xB]);
            C_XX_intensity.D0_only        = C_XX_intensity.D0_only + tmp;
            C_XX_intensity.counts3(:,:,2) = C_XX_intensity.counts3(:,:,2) + tmp;
        end

        maskRest = (oo == 1) | (oo == 4);
        if any(maskRest)
            tmp = accum2_count(ia(maskRest), ib(maskRest), [n_xA, n_xB]);
            C_XX_intensity.rest           = C_XX_intensity.rest + tmp;
            C_XX_intensity.counts3(:,:,3) = C_XX_intensity.counts3(:,:,3) + tmp;
        end
    end
end
    fprintf("本轮发送光子数:%.4f", count);


C_windows = struct('ZZ', C_ZZ, 'XX', C_XX);

end


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

zSend = isZ & (rZ < epsilon);
zVac  = isZ & ~zSend;

mu(zSend)  = signal_mu;
cat(zSend) = 1;

mu(zVac)   = vacuum_mu;
cat(zVac)  = 2;

isX = ~isZ;

% ================= 唯一修改点（GPU安全 + 更快） =================
if any(isX)
    rX_sub = rX(isX);

    % 向量化CDF采样（替代for循环 / find / arrayfun）
    cmp = rX_sub.' <= cdf_decoy;
    idx = sum(~cmp, 2) + 1;
    idx = idx.';

    % 防止极小概率越界
    idx(idx > numel(x_set)) = numel(x_set);

    mu(isX)  = x_set(idx);
    cat(isX) = idx + 2;
end
% =============================================================

end


% =========================================================================
function mu_out = apply_intensity_fluctuation_local_gpu(mu_in, fluct_sigma)

if fluct_sigma <= 0
    mu_out = mu_in;
    return;
end

mu_out = mu_in;
mask = mu_in > 0;

mu_valid = mu_in(mask);

% ---------- 定义区间 ----------
lower = (1 - fluct_sigma) .* mu_valid;
upper = (1 + fluct_sigma) .* mu_valid;

% ---------- 高斯参数 ----------
sigma_abs = fluct_sigma .* mu_valid;

% ---------- 初始化 ----------
mu_sample = mu_valid;

% ---------- GPU / CPU ----------
isGPU = isa(mu_in, 'gpuArray');

% ---------- 拒绝采样 ----------
max_iter = 10;  % 通常2~3次就够

for iter = 1:max_iter
    
    if isGPU
        noise = randn(size(mu_valid), 'like', mu_valid);
    else
        noise = randn(size(mu_valid));
    end

    proposal = mu_valid + sigma_abs .* noise;

    valid = (proposal >= lower) & (proposal <= upper);

    % 接受
    mu_sample(valid) = proposal(valid);

    % 如果全部满足就提前结束
    if all(valid)
        break;
    end

    % 只对不满足的继续采样
    mu_valid = mu_sample;  % 当前值更新（可选）
end

% ---------- 赋值 ----------
mu_out(mask) = mu_sample;

end


% =========================================================================
function c = count_outcomes4(outcome_vec)
if isempty(outcome_vec)
    c = zeros(1, 4);
else
    c = accumarray(double(outcome_vec(:)), 1, [4, 1]).';
end
end


% =========================================================================
function A = accum2_count(i1, i2, out_sz)
if isempty(i1)
    A = zeros(out_sz);
else
    A = accumarray([i1(:), i2(:)], 1, out_sz);
end
end