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
varR = {'signal_A', 'decoy_A', 'signal_B', 'decoy_B'};
varF = setdiff(variables, varR);   % fixed variables

% Pulse and optimization parameters
N = 3.5*1e8;
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
p_decoy = [0.2, 0.6, 0.2];

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
dAC = 50;

signal_A = us;
decoy_A = [u_d1, u_d2];

thetaA = {signal_A, decoy_A, alpha_A, dAC};

% Bob's Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_B = 0.21;
dBC = 50;

signal_B = us;
decoy_B = [u_d1, u_d2];

thetaB = {signal_B, decoy_B, alpha_B, dBC};

% Construct the true parameter set
thetas = {thetaA, thetaB, thetaC};

colSizes = [numel(decoy_A), numel(decoy_B)];

% Compute prior parameters thetaP, and extract true random variables thetaR and fixed variables thetaF
[thetaP, thetaR, thetaF] = processParams(thetas, varF, varR, fluct_sigma);

% -------------------------------------------------------------------------
% Run simulation
% -------------------------------------------------------------------------
    [C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
        simulate_SNS_TF_QKD_3(N, thetaA, thetaB, thetaC, ...
        'buffer', buffer, ...
        'z_window', z_window, ...
        'epsilon', epsilon, ...
        'p_decoy', p_decoy, ...
        'lambda', lambda, ...
        'fluct_sigma', fluct_sigma, ...
        'vacuum', vacuum); 
    disp(C_XX_intensity.counts3);

                % Simulation event order:
            % [D1_only, D0_only, rest]
            Ck = squeeze(C_XX_intensity.counts3(2,2,:)).';
disp(C_XX_intensity.counts3(2,2,:))
            % Compress to:
            % [only D0, only D1, rest]
            C_counts = [Ck(2), Ck(1), Ck(3)];
            % disp(C_counts);