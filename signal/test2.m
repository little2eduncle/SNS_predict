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
             'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1', ...
             'signal_B', 'decoy_B', 'alpha_B', 'dBC'};

% List of random variables
varR_XX = {'decoy_A', 'decoy_B'};
varF_XX = setdiff(variables, varR_XX);   % XX windows fixed variables

% Pulse and optimization parameters
N = 1e7;
algorithm = 'quasi-newton';   % trust-region or quasi-newton
maxIters = 1000;
buffer = 1e4;

% Sampling options
Ns = 1e3;
Nb = 100;
display = false;
chunkSize = 10;
method = 'srss';

% Protocol parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_window = 0.5;
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

% Construct the true parameter set
thetas = {thetaA, thetaB, thetaC};

colSizes = [numel(decoy_A), numel(decoy_B)];

% Compute prior parameters thetaP, and extract true random variables thetaR and fixed variables thetaF
[thetaP_XX, thetaR_XX, thetaF_XX] = processParams(thetas, varF_XX, varR_XX, fluct_sigma);

L_values = 0:0.5:500;
num_runs = numel(L_values);

for idx = 1:num_runs
    l = L_values(idx);
    dAC = l;
    dBC = l;

    thetaA = {signal_A, decoy_A, alpha_A, dAC};
    thetaB = {signal_B, decoy_B, alpha_B, dBC};

    tic

    [C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
        simulate_SNS_TF_QKD_4(N, thetaA, thetaB, thetaC, ...
        'buffer', buffer, ...
        'z_window', z_window, ...
        'epsilon', epsilon, ...
        'p_decoy', p_decoy, ...
        'lambda', lambda, ...
        'fluct_sigma', fluct_sigma, ...
        'vacuum', vacuum);


    [thetaR_MAP, info] = MAP(C_XX_intensity, varR_XX, thetaF_XX, thetaP_XX, colSizes, ...
        'lambda', lambda, ...
        'fluct_sigma', fluct_sigma, ...
        'algorithm', algorithm, ...
        'maxIters', maxIters, ...
        'nStarts', 10, ...
        'relNoise', 0.05, ...
        'doDiagnostics', true);

    t = toc;
    fprintf("第%d轮仿真运行%.1f秒", dAC, t);
    disp(thetaR_MAP);
end