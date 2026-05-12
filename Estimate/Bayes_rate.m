S_ZZ = load('thetaR_ZZ_MAP_results.mat');
S_XX = load('thetaR_XX_MAP_results.mat');

N_m = 1e7;
frac = 0.043;
pc = 0.6;
pd = 1e-8;
alpha = 0.21;
e_det = 0.02;
f = 1.15;
sig = 0.02;
epsilon = 0.275;

R_Bayes = zeros(1, length(S_ZZ.L_values));

for i = 1:length(S_ZZ.L_values)

    d = S_ZZ.L_values(i);

    theta_ZZ = S_ZZ.thetaR_MAP_mat(i,:);
    theta_XX = S_XX.thetaR_MAP_mat(i,:);

    signal_A = theta_ZZ(1);
    signal_B = theta_ZZ(2);

    decoy_A = [theta_XX(1), theta_XX(2)];
    decoy_B = [theta_XX(3), theta_XX(4)];

    [s1_Bayes, e1_Bayes, P] = Z1_data_Bayes(signal_A, decoy_A, signal_B, decoy_B, pd, pc, alpha, d, e_det);

    eta_trans = 10^(-alpha*d/10)*pc;
    S_s = 1 - (1-pd)^2*exp(-2*0.35*eta_trans);
    S_s_half = 1 - (1-pd)^2*exp(-0.35*eta_trans);
    S_0 = 1 - (1-pd)^2;

    S_z = epsilon^2*S_s + (1-epsilon)^2*S_0 + 2*epsilon*(1-epsilon)*S_s_half;
    E_z = (epsilon^2*S_s + (1-epsilon)^2*S_0) / S_z;

    P_s_1_L_Bayes = P{1}(1);
    R_Bayes(i) = max(2*epsilon*(1-epsilon)*P_s_1_L_Bayes*s1_Bayes*(1-h(e1_Bayes))-f*S_z*E_z, 0);
end

% 绘图
figure;
semilogy(S_ZZ.L_values, R_Bayes, 'LineWidth', 1.5);
xlabel('Distance (km)');
ylabel('Secure key rate per pulse');
title('R_{UB}');
grid on;