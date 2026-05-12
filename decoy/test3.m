clear

rng(1);  % optional: fix random seed
addpath(genpath(fileparts(mfilename('fullpath'))));

% Plotting Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale = 2;                    % Factor scaling for screen readability
FigureWidth = 180;            % Width of the figure in mm
FontSize = 8;                 % Font size for plot labels and text
FontName = 'Times New Roman'; % Font name for plot labels and text
Interpreter = 'latex';        % Font rendering (latex or tex)
CIp = 0.99;                   % Confidence interval threshold %#ok<NASGU>

% Scale the figure size and font
FigureWidth = FigureWidth * scale;
FontSize = FontSize * scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define all variable names (order must match unpacking in sampleParameters/processParams)
variables = {'signal_A', 'decoy_A', 'alpha_A', 'dAC', ...
             'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', 'phi_center', ...
             'signal_B', 'decoy_B', 'alpha_B', 'dBC'};

% List of random variables
varR = {'decoy_A', 'decoy_B'};
varF = setdiff(variables, varR);   % fixed variables

% Pulse and optimization parameters
N = 3.9*1e7;
algorithm = 'quasi-newton';   % trust-region or quasi-newton
maxIters = 1000;
buffer = 1e6;

% Sampling options
Ns = 1e3;
Nb = 100;
display = false;
chunkSize = 10;
method = 'srss'; %#ok<NASGU>

% Protocol parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_window = 0.8;
epsilon = 0.978;
p_decoy = [0.6, 0.2, 0.2];

lambda = 0.02;
fluct_sigma = 0.01;

us = 0.35;
u_d1 = 0.01;
u_d2 = 0.12;
vacuum = 0;

% Charlie's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pa = 0;
pc = 0.6;
pd = 1e-7;
phi_center = 0;

% Apply slight variation
pa0 = pa * 0.9;   pa1 = pa * 1.1;
pc0 = pc * 1;     pc1 = pc * 1;
pd0 = pd * 1;     pd1 = pd * 1;

% Group Charlie's parameters
thetaC = {pa0, pa1, pc0, pc1, pd0, pd1, phi_center};

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

% Construct the true parameter set
thetas = {thetaA, thetaB, thetaC};

colSizes = [numel(decoy_A), numel(decoy_B)];

% Compute prior parameters thetaP, and extract true random variables thetaR and fixed variables thetaF
[thetaP, thetaR, thetaF] = processParams(thetas, varF, varR, fluct_sigma);

% -------------------------------------------------------------------------
% Prepare storage for thetaR_MAP
% -------------------------------------------------------------------------
L_values = 0:10:100;
num_runs = numel(L_values);

% thetaR_MAP 预计长度 = sum(colSizes) = 4
map_dim = sum(colSizes);

save_file = 'thetaR_MAP_results.mat';

% 如果文件已存在，先删除，保证每次重新执行都从空文件开始
if exist(save_file, 'file')
    delete(save_file);
end

% 创建 mat 文件，并预分配
M = matfile('thetaR_MAP_results.mat', 'Writable', true);
M.thetaR_MAP_all = nan(num_runs, map_dim);
M.L_values = L_values;
M.colSizes = colSizes;

% 如果你也想保存每次 MAP 的目标函数值/状态，可以一起预分配
M.fval_all = nan(num_runs, 1);

% -------------------------------------------------------------------------
% Run simulation
% -------------------------------------------------------------------------
tic;

for idx = 1:num_runs
    l = L_values(idx);
    dAC = l;
    dBC = l;

    thetaA = {signal_A, decoy_A, alpha_A, dAC};
    thetaB = {signal_B, decoy_B, alpha_B, dBC};

    [C_total, C_windows, C_intensity, C_XX_intensity] = ...
        simulate_SNS_TF_QKD_2(N, thetaA, thetaB, thetaC, ...
        'buffer', buffer, ...
        'z_window', z_window, ...
        'epsilon', epsilon, ...
        'p_decoy', p_decoy, ...
        'lambda', lambda, ...
        'fluct_sigma', fluct_sigma, ...
        'vacuum', vacuum); 

    [thetaR_MAP, info] = MAP(C_XX_intensity, varR, thetaF, thetaP, colSizes, ...
        'lambda', lambda, ...
        'algorithm', algorithm, ...
        'maxIters', maxIters, ...
        'nStarts', 20, ...
        'relNoise', 0.05, ...
        'doDiagnostics', true);
disp(thetaR_MAP);
    % 保证 thetaR_MAP 以 1 x map_dim 行向量写入
    thetaR_MAP_row = thetaR_MAP(:).';

    if numel(thetaR_MAP_row) ~= map_dim
        error('thetaR_MAP length mismatch: expected %d, got %d.', ...
              map_dim, numel(thetaR_MAP_row));
    end

    % 每次循环追加一行
    M.thetaR_MAP_all(idx, :) = thetaR_MAP_row;

    % 可选：保存对应的目标函数值
    if isstruct(info) && isfield(info, 'fval')
        M.fval_all(idx, 1) = info.fval;
    end

    fprintf('Saved thetaR_MAP for idx = %d / %d, L = %.1f km\n', ...
            idx, num_runs, l);
end

t = toc;
disp(['经过时间：', num2str(t), '秒']);