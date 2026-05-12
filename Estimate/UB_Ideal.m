signal_A = 0.496;
signal_B = 0.496;

A1 = 0.1;
A2 = 0.411;
decoy_A = [A1, A2];

B1 = 0.1;
B2 = 0.411;
decoy_B = [B1, B2];

N_m = 1e7;
frac = 0.043;
pc = 0.6;
pd = 1e-8;
alpha = 0.21;
e_det = 0.02;
f = 1.15;
sig = 0.01;
epsilon = 0.275; 

d_list = 0:0.5:500;              % 0 到 500 km，步长 0.5
R_UB = zeros(size(d_list));      % 记录每个距离下的安全码率
R_ideal = zeros(size(d_list));      % 记录每个距离下的安全码率

for i = 1:length(d_list)
    d = d_list(i);

    eta_trans = 10^(-alpha*d/10) * pc;
    S_s = 1 - (1-pd)^2 * exp(-2*0.35*eta_trans);
    S_s_half = 1 - (1-pd)^2 * exp(-0.35*eta_trans);
    S_0 = 1 - (1-pd)^2;

    S_z = epsilon^2 * S_s + (1-epsilon)^2 * S_0 + 2*epsilon*(1-epsilon) * S_s_half;
    E_z = (epsilon^2 * S_s + (1-epsilon)^2 * S_0) / S_z;

    [s1_UB, e1_UB, P_UB] = Z1_data_UB(N_m, frac, sig, pd, pc, alpha, d, e_det);
    [s1_ideal, e1_ideal, P_ideal] = Z1_data_ideal(pd, pc, alpha, d, e_det);

    P_s_1_L_ideal = P_ideal{1}(1);
    P_s_1_L_UB = P_UB{1}(1);

    R_UB(i) = max(2*epsilon*(1-epsilon)*P_s_1_L_UB*s1_UB*(1 - h(e1_UB)) - f*S_z*h(E_z), 0);
    R_ideal(i) = max(2*epsilon*(1-epsilon)*P_s_1_L_ideal*s1_ideal*(1 - h(e1_ideal)) - f*S_z*h(E_z), 0);

end

% ========== 合并绘图 ==========
figure;
semilogy(d_list, R_ideal, 'r-', 'LineWidth', 1.5); hold on;
semilogy(d_list, R_UB, 'b-', 'LineWidth', 1.5);
xlabel('Distance (km)');
ylabel('Secure key rate per pulse');
title('Comparison of R_{UB} and R_{ideal}');
legend('R_{UB}', 'R_{ideal}');
grid on;