clear

data_dir = 'D:\Documents\MATLAB\SNS\sim_data';
save_dir = 'D:\Documents\MATLAB\SNS\predict_data';

% 自动创建保存文件夹
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

rng(1);  % optional: fix random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variables = {'signal_A', 'decoy_A', 'alpha_A', 'dAC', ...
             'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', ...
             'signal_B', 'decoy_B', 'alpha_B', 'dBC'};

varR_XX = {'decoy_A', 'decoy_B'};
varF_XX = setdiff(variables, varR_XX);   % XX windows fixed variables

% Pulse and optimization parameters
N = 3.78*1e7;
algorithm = 'quasi-newton';
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
pc = 0.6;
pd = 1e-8;

pa0 = pa * 0.9;   pa1 = pa * 1.1;
pc0 = pc * 1;     pc1 = pc * 1;
pd0 = pd * 1;     pd1 = pd * 1;

thetaC = {pa0, pa1, pc0, pc1, pd0, pd1, lambda};

% Alice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_A = 0.21;
dAC = 0;

signal_A = us;
decoy_A = [u_d1, u_d2];

thetaA = {signal_A, decoy_A, alpha_A, dAC};

% Bob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_B = 0.21;
dBC = 0;

signal_B = us;
decoy_B = [u_d1, u_d2];

thetaB = {signal_B, decoy_B, alpha_B, dBC};

thetas = {thetaA, thetaB, thetaC};

colSizes = [numel(decoy_A), numel(decoy_B)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== 加载数据 =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_XX = load(fullfile(data_dir, 'C_XX_intensity_all_9.mat'));

L_values = S_XX.L_values;
num_runs = numel(L_values);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====新增：结果存储 =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaR_MAP_all = cell(num_runs, 1);     % 通用存储（推荐）
thetaR_MAP_mat = zeros(num_runs, 4);    % 数值矩阵（decoy_A, decoy_B）

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== 主循环 =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx = 1:num_runs
tic
    L = L_values(idx);
    dAC = L;
    dBC = L;

    thetaA = {signal_A, decoy_A, alpha_A, dAC};
    thetaB = {signal_B, decoy_B, alpha_B, dBC};

    thetas = {thetaA, thetaB, thetaC};

    % 参数处理
    [thetaP_XX, thetaR_XX, thetaF_XX] = ...
    processParams(thetas, varF_XX, varR_XX, fluct_sigma);

    % ---------- XX ----------
    C_XX_counts3 = S_XX.C_XX_counts3_all(:, :, :, idx);

    [thetaR_MAP, info] = MAP_XX(C_XX_counts3, varR_XX, thetaF_XX, thetaP_XX, colSizes, ...
        'lambda',lambda, ...
        'fluct_sigma', fluct_sigma, ...
        'algorithm', algorithm, ...
        'maxIters', maxIters, ...
        'nStarts', 10, ...
        'relNoise', 0.01, ...
        'doDiagnostics', true);
    
    fprintf('L = %.1f km\n', L);
    disp(thetaR_MAP);

  t = toc;
  disp(t);
    % =====新增：记录 =====
    thetaR_MAP_all{idx} = thetaR_MAP;
    thetaR_MAP_mat(idx, :) = thetaR_MAP(:).';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====新增：保存到文件 =====
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(fullfile(save_dir, 'thetaR_XX_MAP_results_7.mat'), ...
    'thetaR_MAP_all', ...
    'thetaR_MAP_mat', ...
    'L_values');

fprintf('结果已保存到 thetaR_XX_MAP_results_7.mat\n');

