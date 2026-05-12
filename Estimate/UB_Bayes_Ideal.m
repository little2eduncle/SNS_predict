% 清空并加载数据
clear; clc;
S_ZZ = load('thetaR_ZZ_MAP_results.mat');
S_XX = load('thetaR_XX_MAP_results.mat');

% ========== 公共参数 ==========
N_m   = 1e7;
frac  = 0.043;
pc    = 0.6;
pd    = 1e-8;          % 统一使用 1e-8（或根据实际需要修改）
alpha = 0.21;
e_det = 0.02;          % 统一使用 0.02
f     = 1.15;
sig   = 0.01;
epsilon = 0.275;

signal_A = 0.496;
signal_B = 0.496;

decoy_A = [0.1, 0.411];
decoy_B = [0.1, 0.411];

d_list = 0:0.5:500;

% ========== 计算 R_ideal ==========
R_ideal   = zeros(size(d_list));

for i = 1:length(d_list)
    d = d_list(i);

    eta_trans = 10^(-alpha*d/10) * pc;
    S_s       = 1 - (1-pd)^2 * exp(-2*signal_A*eta_trans);
    S_s_half  = 1 - (1-pd)^2 * exp(-signal_A*eta_trans);
    S_0       = 1 - (1-pd)^2;

    S_z = epsilon^2 * S_s + (1-epsilon)^2 * S_0 + 2*epsilon*(1-epsilon)*S_s_half;
    E_z = (epsilon^2 * S_s + (1-epsilon)^2 * S_0) / S_z;

    [s1, e1, P] = Z1_data_ideal(signal_A, decoy_A, signal_B, decoy_B, pd, pc, alpha, d, e_det);
    P_s_1_L = P{1}(1);

    R_ideal(i)  = max(2*epsilon*(1-epsilon)*P_s_1_L*s1*(1 - h(e1)) - f*S_z*h(E_z), 0);
end

% ========== 计算 R_UB ==========
R_UB   = zeros(size(d_list));
s1_UB  = zeros(size(d_list));
e1_UB  = zeros(size(d_list));

for i = 1:length(d_list)
    d = d_list(i);

    eta_trans = 10^(-alpha*d/10) * pc;
    S_s       = 1 - (1-pd)^2 * exp(-2*signal_A*eta_trans);
    S_s_half  = 1 - (1-pd)^2 * exp(-signal_A*eta_trans);
    S_0       = 1 - (1-pd)^2;

    S_z = epsilon^2 * S_s + (1-epsilon)^2 * S_0 + 2*epsilon*(1-epsilon)*S_s_half;
    E_z = (epsilon^2 * S_s + (1-epsilon)^2 * S_0) / S_z;

    [s1, e1, P] = Z1_data_UB(signal_A, decoy_A, signal_B, decoy_B, N_m, frac, sig, pd, pc, alpha, d, e_det);

    P_s_1_L = P{1}(1);
    R_UB(i)  = max(2*epsilon*(1-epsilon)*P_s_1_L*s1*(1 - h(e1)) - f*S_z*h(E_z), 0);
end

% ========== 计算 R_Bayes ==========
R_Bayes = zeros(1, length(S_ZZ.L_values));

for i = 1:length(S_ZZ.L_values)
    d = S_ZZ.L_values(i);

    theta_ZZ = S_ZZ.thetaR_MAP_mat(i, :);
    theta_XX = S_XX.thetaR_MAP_mat(i, :);

    signal_A_predict = theta_ZZ(1);
    signal_B_predict = theta_ZZ(2);
    decoy_A_predict  = [theta_XX(1), theta_XX(2)];
    decoy_B_predict  = [theta_XX(3), theta_XX(4)];

    [s1_Bayes, e1_Bayes, P] = Z1_data_Bayes2(...
        signal_A_predict, decoy_A_predict, signal_B_predict, decoy_B_predict, pd, pc, alpha, d, e_det);

    eta_trans = 10^(-alpha*d/10)*pc;

    S_0 = 1 - (1-pd)^2;

    S_ua_0 = 1 - (1-pd)^2*exp(-signal_A_predict*eta_trans);
    S_0_ub = 1 - (1-pd)^2*exp(-signal_B_predict*eta_trans);
    S_s_half = S_ua_0 + S_0_ub;

    S_s = 1 - (1-pd)^2*exp(-(signal_A_predict+signal_B_predict)*eta_trans);

    S_z = epsilon^2*S_s + (1-epsilon)^2*S_0 + epsilon*(1-epsilon)*S_s_half;
    E_z = (epsilon^2*S_s + (1-epsilon)^2*S_0) / S_z;


    P_s_1_L_Bayes = P{1}(1);

    R_Bayes(i) = max(2*epsilon*(1-epsilon)*P_s_1_L_Bayes*s1_Bayes*(1 - h(e1_Bayes)) - f*S_z*h(E_z), 0);
end

% ========== 合并绘图 ==========
figure;
semilogy(d_list, R_ideal, 'k-', 'LineWidth', 1.5); hold on;
semilogy(d_list, R_UB, 'b-', 'LineWidth', 1.5); hold on;
semilogy(S_ZZ.L_values, R_Bayes, 'r--', 'LineWidth', 1.5);
xlabel('Distance (km)');
ylabel('Secure key rate per pulse');
title('Comparison of R_{ideal} ,R_{UB} and R_{Bayes}');
legend({'R_{ideal}', 'R_{UB}', 'R_{Bayes}'}, 'Location', 'northeast');
grid on;