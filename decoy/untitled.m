clear

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
N = 1e7;
algorithm = 'quasi-newton';
maxIters = 1000;
buffer = 1e4;

% Protocol parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_window = 0.866;
epsilon = 0.275;
p_decoy = [0.898, 0.066, 0.036];

lambda = 0.015;
fluct_sigma = 0.01;

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
dAC = 500;

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
%==============================================
colSizes = [numel(decoy_A), numel(decoy_B)];

[thetaP_XX, thetaR_XX, thetaF_XX] = ...
    processParams(thetas, varF_XX, varR_XX, fluct_sigma);

[C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
    simulate_SNS_TF_QKD_5(N, thetaA, thetaB, thetaC, ...
    'buffer', buffer, ...
    'z_window', z_window, ...
    'epsilon', epsilon, ...
    'p_decoy', p_decoy, ...
    'lambda', lambda, ...
    'fluct_sigma', fluct_sigma, ...
    'vacuum', vacuum);

disp(C_XX_intensity.counts3);

[thetaR_MAP, info] = MAP(C_XX_intensity.counts3, varR_XX, thetaF_XX, thetaP_XX, colSizes, ...
    'fluct_sigma', fluct_sigma, ...
    'algorithm', algorithm, ...
    'maxIters', maxIters, ...
    'nStarts', 10, ...
    'relNoise', 0.01, ...
    'doDiagnostics', true);

disp(thetaR_MAP);

