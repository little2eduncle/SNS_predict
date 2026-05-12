clear

% rng(1);  % optional: fix random seed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Session Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

variables = {'signal_A', 'decoy_A', 'alpha_A', 'dAC', ...
             'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', ...
             'signal_B', 'decoy_B', 'alpha_B', 'dBC'};

varR_ZZ = {'signal_A', 'signal_B'};
varF_ZZ = setdiff(variables, varR_ZZ);

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
dAC = 200;

signal_A = us;
decoy_A = [u_d1, u_d2];

thetaA = {signal_A, decoy_A, alpha_A, dAC};

% Bob
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_B = 0.21;
dBC = 200;

signal_B = us;
decoy_B = [u_d1, u_d2];

thetaB = {signal_B, decoy_B, alpha_B, dBC};

thetas = {thetaA, thetaB, thetaC};
%==============================================
colSizes = [1, 1];

[thetaP_ZZ, thetaR_ZZ, thetaF_ZZ] = ...
    processParams(thetas, varF_ZZ, varR_ZZ, fluct_sigma);

tic;

[C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
    simulate_SNS_TF_QKD_5(N, thetaA, thetaB, thetaC, ...
    'buffer', buffer, ...
    'z_window', z_window, ...
    'epsilon', epsilon, ...
    'p_decoy', p_decoy, ...
    'lambda', lambda, ...
    'fluct_sigma', fluct_sigma, ...
    'vacuum', vacuum);

disp(C_ZZ_intensity.counts3);
disp("-----");
disp(C_XX_intensity.counts3);

t = toc;
disp(t);


