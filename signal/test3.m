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

varR_ZZ = {'signal_A', 'signal_B'};
varF_ZZ = setdiff(variables, varR_ZZ);   % ZZ windows fixed variables

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

colSizes = [1, 1];

S_ZZ = load('C_ZZ_intensity_all.mat');

L_values = S_ZZ.L_values;
num_runs = numel(L_values);

for idx = 1:num_runs
    L = L_values(idx);
    dAC = L;
    dBC = L;
    thetaA = {signal_A, decoy_A, alpha_A, dAC};
    thetaB = {signal_B, decoy_B, alpha_B, dBC};

    thetas = {thetaA, thetaB, thetaC};

    % Compute prior parameters thetaP, and extract true random variables thetaR and fixed variables thetaF
    [thetaP_ZZ, thetaR_ZZ, thetaF_ZZ] = processParams(thetas, varF_ZZ, varR_ZZ, fluct_sigma);

    % ---------- ZZ ----------
    C_ZZ_counts3 = S_ZZ.C_ZZ_counts3_all(:, :, :, idx);

    fprintf('L = %.1f km\n', L);

    [thetaR_MAP, info] = MAP(C_ZZ_counts3, varR_ZZ, thetaF_ZZ, thetaP_ZZ, colSizes, ...
        'fluct_sigma', fluct_sigma, ...
        'algorithm', algorithm, ...
        'maxIters', maxIters, ...
        'nStarts', 10, ...
        'relNoise', 0.01, ...
        'doDiagnostics', true);
    disp(thetaR_MAP);

end
