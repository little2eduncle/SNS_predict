function [C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
    simulate_SNS_TF_QKD_3(N, thetaA, thetaB, thetaC, varargin)
% simulate_SNS_TF_QKD_3
%
% 功能:
%   1) 记录 ZZ 窗口下“所有”光强组合的三类事件 {01, 10, rest}
%   2) 记录 XX 窗口下“通过相位筛选”的所有光强组合的三类事件 {01, 10, rest}
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
%   'fluct_sigma' : 非零光强的高斯波动标准差
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
%       Alice 的 Z 态:
%           1 -> signal_A
%           2 -> vacuum
%       Bob 的 Z 态:
%           1 -> signal_B
%           2 -> vacuum
%
%       字段:
%           .mu_nominal_A   : [2x1]，Alice 的名义光强 [signal_A; vacuum]
%           .mu_nominal_B   : [1x2]，Bob 的名义光强 [signal_B, vacuum]
%           .counts3        : [2x2x3]
%                             counts3(i,j,1) = 01 次数
%                             counts3(i,j,2) = 10 次数
%                             counts3(i,j,3) = rest={00,11} 次数
%           .D1_only        : [2x2]，对应 counts3(:,:,1)
%           .D0_only        : [2x2]，对应 counts3(:,:,2)
%           .rest           : [2x2]，对应 counts3(:,:,3)
%           .total_valid    : [2x2]，该 ZZ 光强组合总出现次数
%
%   C_XX_intensity : struct
%       XX 窗口下所有 X 光强组合的三类事件统计，
%       但只统计“通过相位筛选”的事件
%
%       字段:
%           .mu_nominal_A   : [n_xA x 1]，Alice 的 X 窗口名义光强
%           .mu_nominal_B   : [1 x n_xB]，Bob 的 X 窗口名义光强
%           .counts3        : [n_xA x n_xB x 3]
%                             counts3(i,j,1) = 01 次数
%                             counts3(i,j,2) = 10 次数
%                             counts3(i,j,3) = rest={00,11} 次数
%           .D1_only        : [n_xA x n_xB]
%           .D0_only        : [n_xA x n_xB]
%           .rest           : [n_xA x n_xB]
%           .total_valid    : [n_xA x n_xB]
%                             表示该 XX 光强组合中，通过相位筛选后的总次数

% -------------------------------------------------------------------------
% 解析输入参数
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'buffer', 1e4);        % 每个 chunk 处理的样本数
addParameter(p, 'z_window', 0.5);      % 选择 Z 窗口的概率
addParameter(p, 'epsilon', 0.99);      % Z 窗口中发送 signal 的概率
addParameter(p, 'p_decoy', []);        % X 窗口下各 decoy/vacuum 的概率
addParameter(p, 'lambda', 0.2);        % 相位筛选参数
addParameter(p, 'fluct_sigma', 0.0);   % 光强高斯波动标准差
addParameter(p, 'vacuum', 0);          % 真空态光强
addParameter(p, 'use_gpu', true);      % 是否使用 GPU
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

% 如果要求使用 GPU，则检查 GPU 是否可用
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

% 保证 decoy 是行向量
decoy_A = decoy_A(:).';
decoy_B = decoy_B(:).';

% X 窗口可选光强集合 = [decoys..., vacuum]
x_set_A = [decoy_A, vacuum];
x_set_B = [decoy_B, vacuum];

% 如果没有指定 p_decoy，则对 X 窗口各个光强均匀选取
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
n_signalA = 2;   % Alice 的 Z 态数：signal / vacuum
n_signalB = 2;   % Bob   的 Z 态数：signal / vacuum
n_xA      = numel(x_set_A);
n_xB      = numel(x_set_B);

% 所有脉冲的 outcome 统计 [00, 01, 10, 11]
C_total = zeros(1, 4);

% 窗口级别统计
C_ZZ = zeros(1, 4);
C_XX = zeros(1, 4);

% -----------------------------
% ZZ 窗口下所有 2x2 光强组合
% -----------------------------
C_ZZ_intensity.mu_nominal_A = [signal_A; vacuum];
C_ZZ_intensity.mu_nominal_B = [signal_B, vacuum];

% 三类计数:
%   1 -> 01
%   2 -> 10
%   3 -> rest = {00,11}
C_ZZ_intensity.counts3     = zeros(n_signalA, n_signalB, 3);

% 单独存这三类，便于后续直接调用
C_ZZ_intensity.D1_only     = zeros(n_signalA, n_signalB);
C_ZZ_intensity.D0_only     = zeros(n_signalA, n_signalB);
C_ZZ_intensity.rest        = zeros(n_signalA, n_signalB);

% 每个 ZZ 组合出现的总次数
C_ZZ_intensity.total_valid = zeros(n_signalA, n_signalB);

% ---------------------------------------
% XX 窗口下所有光强组合（仅相位筛选通过）
% ---------------------------------------
C_XX_intensity.mu_nominal_A = x_set_A(:);
C_XX_intensity.mu_nominal_B = x_set_B(:).';

% 三类计数:
%   1 -> 01
%   2 -> 10
%   3 -> rest = {00,11}
C_XX_intensity.counts3     = zeros(n_xA, n_xB, 3);

C_XX_intensity.D1_only     = zeros(n_xA, n_xB);
C_XX_intensity.D0_only     = zeros(n_xA, n_xB);
C_XX_intensity.rest        = zeros(n_xA, n_xB);

% 每个 XX 组合在“通过相位筛选后”的总次数
C_XX_intensity.total_valid = zeros(n_xA, n_xB);

% chunk 数
num_chunks = ceil(N / buffer);

% -------------------------------------------------------------------------
% 固定信道透过率
% -------------------------------------------------------------------------
etaA = 10^(-alpha_A * dAC / 10);
etaB = 10^(-alpha_B * dBC / 10);

% 将常数搬到 GPU
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

% 用户如果有 start_up 脚本，这里执行
if exist('start_up', 'file') == 2
    start_up;
end

% -------------------------------------------------------------------------
% 主循环：分块仿真
% -------------------------------------------------------------------------
for chunk = 1:num_chunks
    first = (chunk - 1) * buffer + 1;
    last  = min(chunk * buffer, N);
    chunk_size = last - first + 1;

    % -------------------------------------------------------------
    % 1) 为 Alice / Bob 随机选择窗口：Z 或 X
    % -------------------------------------------------------------
    if use_gpu
        A_isZ = rand(1, chunk_size, 'gpuArray') < z_window;
        B_isZ = rand(1, chunk_size, 'gpuArray') < z_window;
    else
        A_isZ = rand(1, chunk_size) < z_window;
        B_isZ = rand(1, chunk_size) < z_window;
    end

    % -------------------------------------------------------------
    % 2) 根据窗口选择具体发射态，并编码为 catA / catB
    %    Z 窗口:
    %       1 = signal
    %       2 = vacuum
    %    X 窗口:
    %       3..(2+n_x) = x_set 中的光强索引
    % -------------------------------------------------------------
    [muA_nominal, catA] = prepare_party_states_cat( ...
        A_isZ, signal_A, vacuum, x_set_A_dev, epsilon, cdf_decoy_dev, use_gpu);

    [muB_nominal, catB] = prepare_party_states_cat( ...
        B_isZ, signal_B, vacuum, x_set_B_dev, epsilon, cdf_decoy_dev, use_gpu);

    % -------------------------------------------------------------
    % 3) 随机生成相位差 delta_theta ~ Uniform(0, 2*pi)
    % -------------------------------------------------------------
    if use_gpu
        delta_theta = 2 * pi * rand(1, chunk_size, 'gpuArray');
    else
        delta_theta = 2 * pi * rand(1, chunk_size);
    end

    % -------------------------------------------------------------
    % 4) 对非零光强加入高斯波动
    % -------------------------------------------------------------
    muA = apply_intensity_fluctuation_local_gpu(muA_nominal, fluct_sigma);
    muB = apply_intensity_fluctuation_local_gpu(muB_nominal, fluct_sigma);

    % -------------------------------------------------------------
    % 5) 经过信道衰减后的有效光强
    % -------------------------------------------------------------
    muA_eff = muA .* etaA;
    muB_eff = muB .* etaB;

    % -------------------------------------------------------------
    % 6) 在 Charlie 处干涉，得到 D0 / D1 的输入强度 I0 / I1
    % -------------------------------------------------------------
    cdt = cos(delta_theta);
    inter = 2 .* sqrt(muA_eff .* muB_eff) .* cdt;

    I0 = 0.5 .* (muA_eff + muB_eff + inter);
    I1 = 0.5 .* (muA_eff + muB_eff - inter);

    % 数值保护，防止负数
    I0 = max(I0, 0);
    I1 = max(I1, 0);

    % -------------------------------------------------------------
    % 7) 根据暗计数 / 后脉冲 / 探测系数，计算 D0 / D1 点击概率
    % -------------------------------------------------------------
    p0_final = 1 - (1 - pd0) .* (1 - pa0) .* exp(-pc0 .* I0);
    p1_final = 1 - (1 - pd1) .* (1 - pa1) .* exp(-pc1 .* I1);

    % 数值裁剪到 [0,1]
    p0_final = min(max(p0_final, 0), 1);
    p1_final = min(max(p1_final, 0), 1);

    % -------------------------------------------------------------
    % 8) 按点击概率进行伯努利抽样，得到本 chunk 的点击结果
    % -------------------------------------------------------------
    if use_gpu
        finalD0_chunk = rand(1, chunk_size, 'gpuArray') < p0_final;
        finalD1_chunk = rand(1, chunk_size, 'gpuArray') < p1_final;
    else
        finalD0_chunk = rand(1, chunk_size) < p0_final;
        finalD1_chunk = rand(1, chunk_size) < p1_final;
    end

    % outcome 编码:
    %   1 -> 00
    %   2 -> 01 (only D1)
    %   3 -> 10 (only D0)
    %   4 -> 11
    outcome = 1 + double(finalD1_chunk) + 2 * double(finalD0_chunk);

    % -------------------------------------------------------------
    % 9) 定义 ZZ / XX 统计口径
    %
    % ZZ:
    %   所有 catA<=2 且 catB<=2 的组合都记
    %
    % XX:
    %   只记录 catA>=3 且 catB>=3 且通过相位筛选的事件
    % -------------------------------------------------------------
    phase_ok = abs(cdt) >= (1 - lambda);
    valid_xx = (catA >= 3) & (catB >= 3) & phase_ok;
    valid_zz = (catA <= 2) & (catB <= 2);

    % -------------------------------------------------------------
    % 10) 只 gather 小矩阵，避免把大数组全部搬回 CPU
    % -------------------------------------------------------------
    packed = gather(uint8([catA; catB; outcome; valid_zz; valid_xx]));

    A   = packed(1, :);        % Alice 类别编码
    B   = packed(2, :);        % Bob   类别编码
    O   = packed(3, :);        % outcome
    VZZ = packed(4, :) > 0;    % 是否计入 ZZ
    VXX = packed(5, :) > 0;    % 是否计入 XX

    % -------------------------------------------------------------
    % 11) 更新总体 / 窗口级别 outcome 统计
    % -------------------------------------------------------------
    C_total = C_total + count_outcomes4(O);
    C_ZZ    = C_ZZ + count_outcomes4(O(VZZ));
    C_XX    = C_XX + count_outcomes4(O(VXX));

    % -------------------------------------------------------------
    % 12) 统计 ZZ 窗口下各光强组合的三类事件 {01,10,rest}
    % -------------------------------------------------------------
    if any(VZZ)
        ia = double(A(VZZ));    % 1..2
        ib = double(B(VZZ));    % 1..2
        oo = double(O(VZZ));

        % 该 ZZ 组合的总出现次数
        C_ZZ_intensity.total_valid = C_ZZ_intensity.total_valid + ...
            accum2_count(ia, ib, [n_signalA, n_signalB]);

        % 01: only D1
        maskD1 = (oo == 2);
        if any(maskD1)
            tmp = accum2_count(ia(maskD1), ib(maskD1), [n_signalA, n_signalB]);
            C_ZZ_intensity.D1_only        = C_ZZ_intensity.D1_only + tmp;
            C_ZZ_intensity.counts3(:,:,1) = C_ZZ_intensity.counts3(:,:,1) + tmp;
        end

        % 10: only D0
        maskD0 = (oo == 3);
        if any(maskD0)
            tmp = accum2_count(ia(maskD0), ib(maskD0), [n_signalA, n_signalB]);
            C_ZZ_intensity.D0_only        = C_ZZ_intensity.D0_only + tmp;
            C_ZZ_intensity.counts3(:,:,2) = C_ZZ_intensity.counts3(:,:,2) + tmp;
        end

        % rest = {00,11}
        maskRest = (oo == 1) | (oo == 4);
        if any(maskRest)
            tmp = accum2_count(ia(maskRest), ib(maskRest), [n_signalA, n_signalB]);
            C_ZZ_intensity.rest           = C_ZZ_intensity.rest + tmp;
            C_ZZ_intensity.counts3(:,:,3) = C_ZZ_intensity.counts3(:,:,3) + tmp;
        end
    end

    % -------------------------------------------------------------
    % 13) 统计 XX 窗口下各光强组合的三类事件 {01,10,rest}
    %     注意这里只统计“通过相位筛选”的 XX 事件
    % -------------------------------------------------------------
    if any(VXX)
        ia = double(A(VXX) - 2);   % 1..n_xA
        ib = double(B(VXX) - 2);   % 1..n_xB
        oo = double(O(VXX));

        % 该 XX 组合在通过相位筛选后的总出现次数
        C_XX_intensity.total_valid = C_XX_intensity.total_valid + ...
            accum2_count(ia, ib, [n_xA, n_xB]);

        % 01: only D1
        maskD1 = (oo == 2);
        if any(maskD1)
            tmp = accum2_count(ia(maskD1), ib(maskD1), [n_xA, n_xB]);
            C_XX_intensity.D1_only        = C_XX_intensity.D1_only + tmp;
            C_XX_intensity.counts3(:,:,1) = C_XX_intensity.counts3(:,:,1) + tmp;
        end

        % 10: only D0
        maskD0 = (oo == 3);
        if any(maskD0)
            tmp = accum2_count(ia(maskD0), ib(maskD0), [n_xA, n_xB]);
            C_XX_intensity.D0_only        = C_XX_intensity.D0_only + tmp;
            C_XX_intensity.counts3(:,:,2) = C_XX_intensity.counts3(:,:,2) + tmp;
        end

        % rest = {00,11}
        maskRest = (oo == 1) | (oo == 4);
        if any(maskRest)
            tmp = accum2_count(ia(maskRest), ib(maskRest), [n_xA, n_xB]);
            C_XX_intensity.rest           = C_XX_intensity.rest + tmp;
            C_XX_intensity.counts3(:,:,3) = C_XX_intensity.counts3(:,:,3) + tmp;
        end
    end
end

% -------------------------------------------------------------------------
% 窗口级别输出
% -------------------------------------------------------------------------
C_windows = struct('ZZ', C_ZZ, ...
                   'XX', C_XX);

end

% =========================================================================
% prepare_party_states_cat
%
% 根据是否选择 Z 窗口，为每个脉冲生成：
%   1) 名义光强 mu
%   2) 类别编码 cat
%
% 编码规则:
%   Z 窗口:
%       1 = signal
%       2 = vacuum
%   X 窗口:
%       3..(2+n_x) = x_set 中的索引
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

% -----------------------------
% Z 窗口:
%   以 epsilon 概率发送 signal
%   以 1-epsilon 概率发送 vacuum
% -----------------------------
zSend = isZ & (rZ < epsilon);
zVac  = isZ & ~zSend;

mu(zSend)  = signal_mu;
cat(zSend) = 1;

mu(zVac)   = vacuum_mu;
cat(zVac)  = 2;

% -----------------------------
% X 窗口:
%   按 p_decoy 从 x_set = [decoys..., vacuum] 中选取
% -----------------------------
isX = ~isZ;

lower = 0;
for k = 1:(numel(x_set) - 1)
    mask = isX & (rX > lower) & (rX <= cdf_decoy(k));
    mu(mask)  = x_set(k);
    cat(mask) = k + 2;
    lower = cdf_decoy(k);
end

% 最后一档直接接收剩余概率
maskLast = isX & (cat == 0);
mu(maskLast)  = x_set(end);
cat(maskLast) = numel(x_set) + 2;

end

% =========================================================================
% apply_intensity_fluctuation_local_gpu
%
% 对非零光强加入高斯波动:
%   mu_out = mu_in + N(0, fluct_sigma^2)
% 然后裁剪到 >= 0
%
% 真空态 (mu=0) 不加波动
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
% count_outcomes4
%
% 对 outcome = 1/2/3/4 进行计数，返回:
%   [N_00, N_01, N_10, N_11]
% =========================================================================
function c = count_outcomes4(outcome_vec)
if isempty(outcome_vec)
    c = zeros(1, 4);
else
    c = accumarray(double(outcome_vec(:)), 1, [4, 1]).';
end
end

% =========================================================================
% accum2_count
%
% 二维组合计数:
%   A(i1, i2) = 该组合出现次数
%
% 例如:
%   i1 = Alice 的光强索引
%   i2 = Bob   的光强索引
% =========================================================================
function A = accum2_count(i1, i2, out_sz)
if isempty(i1)
    A = zeros(out_sz);
else
    A = accumarray([i1(:), i2(:)], 1, out_sz);
end
end