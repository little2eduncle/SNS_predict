clear

save_dir = 'D:\Documents\MATLAB\SNS\sim_data';

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

rng(1);  % optional: fix random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define all variable names (order must match unpacking in sampleParameters/processParams)
variables = {'signal_A', 'decoy_A', 'alpha_A', 'dAC', ...
             'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', ...
             'signal_B', 'decoy_B', 'alpha_B', 'dBC'};

% Pulse and optimization parameters
N = 3.9*1e7;%3.78*1e7
algorithm = 'quasi-newton';   % trust-region or quasi-newton
maxIters = 1000;
buffer = 1e6;

% Protocol parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_window = 0.866;
epsilon = 0.275;
p_decoy = [0.898, 0.066, 0.036];

lambda = 0.015;
fluct_sigma = 0.05;

us = 0.496;
u_d1 = 0.1;
u_d2 = 0.411;
vacuum = 0;

% Charlie's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pa = 0;
pc = 0.8;
pd = 1e-8;

% Apply slight variation
pa0 = pa * 0.9;   pa1 = pa * 1.1;
pc0 = pc * 1;   pc1 = pc * 1;
pd0 = pd * 1;   pd1 = pd * 1;

% Group Charlie's parameters
thetaC = {pa0, pa1, pc0, pc1, pd0, pd1, lambda};

% Alice's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_A = 0.21;
dAC = 0;

signal_A = us;
decoy_A = [u_d1, u_d2];

thetaA = {signal_A, decoy_A, alpha_A, dAC};

% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_B = 0.21;
dBC = 0;

signal_B = us;
decoy_B = [u_d1, u_d2];

thetaB = {signal_B, decoy_B, alpha_B, dBC};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colSizes = [numel(decoy_A), numel(decoy_B)];

L_values = 0:0.5:250;
num_runs = numel(L_values);

% =========================
% 新增：预分配存储 ZZ / XX intensity
% =========================
n_signalA = 2;
n_signalB = 2;

C_ZZ_counts3_all     = zeros(n_signalA, n_signalB, 3, num_runs);
C_ZZ_D1_only_all     = zeros(n_signalA, n_signalB, num_runs);
C_ZZ_D0_only_all     = zeros(n_signalA, n_signalB, num_runs);
C_ZZ_rest_all        = zeros(n_signalA, n_signalB, num_runs);
C_ZZ_total_valid_all = zeros(n_signalA, n_signalB, num_runs);

C_ZZ_mu_nominal_A = [signal_A; vacuum];
C_ZZ_mu_nominal_B = [signal_B, vacuum];

n_xA = numel([decoy_A, vacuum]);
n_xB = numel([decoy_B, vacuum]);

C_XX_counts3_all     = zeros(n_xA, n_xB, 3, num_runs);
C_XX_D1_only_all     = zeros(n_xA, n_xB, num_runs);
C_XX_D0_only_all     = zeros(n_xA, n_xB, num_runs);
C_XX_rest_all        = zeros(n_xA, n_xB, num_runs);
C_XX_total_valid_all = zeros(n_xA, n_xB, num_runs);
C_XX_total_raw_all   = zeros(n_xA, n_xB, num_runs);

C_XX_mu_nominal_A = [decoy_A, vacuum].';
C_XX_mu_nominal_B = [decoy_B, vacuum];

tic

for idx = 1:num_runs
    l = L_values(idx);
    dAC = l;
    dBC = l;

    thetaA = {signal_A, decoy_A, alpha_A, dAC};
    thetaB = {signal_B, decoy_B, alpha_B, dBC};

    tic

    [C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
        simulate_SNS_TF_QKD_6(N, thetaA, thetaB, thetaC, ...
        'buffer', buffer, ...
        'z_window', z_window, ...
        'epsilon', epsilon, ...
        'p_decoy', p_decoy, ...
        'lambda', lambda, ...
        'fluct_sigma', fluct_sigma, ...
        'vacuum', vacuum);

    t1 = toc;
    fprintf("第%d轮数据生成时间:%f\n", idx, t1);

    % =========================
    % 新增：存储当前距离下的 ZZ / XX intensity
    % =========================
    C_ZZ_counts3_all(:, :, :, idx)     = C_ZZ_intensity.counts3;
    C_ZZ_D1_only_all(:, :, idx)        = C_ZZ_intensity.D1_only;
    C_ZZ_D0_only_all(:, :, idx)        = C_ZZ_intensity.D0_only;
    C_ZZ_rest_all(:, :, idx)           = C_ZZ_intensity.rest;
    C_ZZ_total_valid_all(:, :, idx)    = C_ZZ_intensity.total_valid;

    C_XX_counts3_all(:, :, :, idx)     = C_XX_intensity.counts3;
    C_XX_D1_only_all(:, :, idx)        = C_XX_intensity.D1_only;
    C_XX_D0_only_all(:, :, idx)        = C_XX_intensity.D0_only;
    C_XX_rest_all(:, :, idx)           = C_XX_intensity.rest;
    C_XX_total_valid_all(:, :, idx)    = C_XX_intensity.total_valid;
    C_XX_total_raw_all(:, :, idx)      = C_XX_intensity.total_raw;

end

% =========================
% 新增：循环结束后保存成两个大文件
% =========================
save(fullfile(save_dir, 'C_ZZ_intensity_all_10.mat'), ...
    'C_ZZ_counts3_all', ...
    'C_ZZ_D1_only_all', 'C_ZZ_D0_only_all', 'C_ZZ_rest_all', ...
    'C_ZZ_total_valid_all', ...
    'C_ZZ_mu_nominal_A', 'C_ZZ_mu_nominal_B', ...
    'L_values', '-v7.3');

save(fullfile(save_dir, 'C_XX_intensity_all_10.mat'), ...
    'C_XX_counts3_all', ...
    'C_XX_D1_only_all', 'C_XX_D0_only_all', 'C_XX_rest_all', ...
    'C_XX_total_valid_all', 'C_XX_total_raw_all', ...
    'C_XX_mu_nominal_A', 'C_XX_mu_nominal_B', ...
    'L_values', '-v7.3');

t = toc;
fprintf("数据生成时间:%f秒", t);