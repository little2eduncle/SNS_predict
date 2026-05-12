clear

% -------------------------------------------------------------------------
% 波动参数（非理想光源）
% -------------------------------------------------------------------------
signal_A = 0.35;
signal_B = 0.35;
alpha    = 0.21;
epsilon  = 0.978;
e_det    = 0.02;
f        = 1.1;
pc       = 0.6;
pd       = 1e-7;

d_values = 0:10:500;
R1_all   = nan(size(d_values));
R_ideal_all = nan(size(d_values));

M = matfile('thetaR_MAP_results.mat');

% =========================================================================
% 1) 不可信光源曲线
% =========================================================================
for idx = 1:numel(d_values)

    theta = M.thetaR_MAP_all(idx, :);
    d = d_values(idx);

    % ---------------------------------------------------------------------
    % 预测出来的 decoy 强度
    % ---------------------------------------------------------------------
    decoy_A1_L = theta(1);
    decoy_A2_L = theta(2);
    decoy_B1_L = theta(3);
    decoy_B2_L = theta(4);

    decoy_A1_U = theta(1);
    decoy_A2_U = theta(2);
    decoy_B1_U = theta(3);
    decoy_B2_U = theta(4);

    % ---------------------------------------------------------------------
    % 泊松项
    % ---------------------------------------------------------------------
    P_us_1 = signal_A * exp(-signal_A) * exp(-signal_B) + ...
             signal_B * exp(-signal_B) * exp(-signal_A);

    P_u_d1_0_L = exp(-decoy_A1_L) * exp(-decoy_B1_L);
    P_u_d2_0_L = exp(-decoy_A2_L) * exp(-decoy_B2_L);

    P_u_d1_0_U = exp(-decoy_A1_U) * exp(-decoy_B1_U);
    P_u_d2_0_U = exp(-decoy_A2_U) * exp(-decoy_B2_U);

    P_u_d1_1_L = decoy_A1_L * exp(-decoy_A1_L) * exp(-decoy_B1_L) + ...
                 decoy_B1_L * exp(-decoy_B1_L) * exp(-decoy_A1_L);

    P_u_d2_1_L = decoy_A2_L * exp(-decoy_A2_L) * exp(-decoy_B2_L) + ...
                 decoy_B2_L * exp(-decoy_B2_L) * exp(-decoy_A2_L);

    P_u_d1_1_U = decoy_A1_U * exp(-decoy_A1_U) * exp(-decoy_B1_U) + ...
                 decoy_B1_U * exp(-decoy_B1_U) * exp(-decoy_A1_U);

    P_u_d2_1_U = decoy_A2_U * exp(-decoy_A2_U) * exp(-decoy_B2_U) + ...
                 decoy_B2_U * exp(-decoy_B2_U) * exp(-decoy_A2_U);

    P_u_d1_2_L = exp(-decoy_A1_L) * 0.5 * decoy_B1_L^2 * exp(-decoy_B1_L) + ...
                 exp(-decoy_B1_L) * 0.5 * decoy_A1_L^2 * exp(-decoy_A1_L) + ...
                 decoy_A1_L * exp(-decoy_A1_L) * decoy_B1_L * exp(-decoy_B1_L);

    P_u_d2_2_L = exp(-decoy_A2_L) * 0.5 * decoy_B2_L^2 * exp(-decoy_B2_L) + ...
                 exp(-decoy_B2_L) * 0.5 * decoy_A2_L^2 * exp(-decoy_A2_L) + ...
                 decoy_A2_L * exp(-decoy_A2_L) * decoy_B2_L * exp(-decoy_B2_L);

    P_u_d1_2_U = exp(-decoy_A1_U) * 0.5 * decoy_B1_U^2 * exp(-decoy_B1_U) + ...
                 exp(-decoy_B1_U) * 0.5 * decoy_A1_U^2 * exp(-decoy_A1_U) + ...
                 decoy_A1_U * exp(-decoy_A1_U) * decoy_B1_U * exp(-decoy_B1_U);

    P_u_d2_2_U = exp(-decoy_A2_U) * 0.5 * decoy_B2_U^2 * exp(-decoy_B2_U) + ...
                 exp(-decoy_B2_U) * 0.5 * decoy_A2_U^2 * exp(-decoy_A2_U) + ...
                 decoy_A2_U * exp(-decoy_A2_U) * decoy_B2_U * exp(-decoy_B2_U);

    % ---------------------------------------------------------------------
    % 增益和误码率
    % ---------------------------------------------------------------------
    eta = 10^(-alpha * d / 10);

    us   = signal_A + signal_B;
    u_d1 = theta(1) + theta(3);
    u_d2 = theta(2) + theta(4);

    S_us   = 1 - (1 - pd)^2 * exp(-us   * pc * eta);
    S_u_d1 = 1 - (1 - pd)^2 * exp(-u_d1 * pc * eta);
    S_u_d2 = 1 - (1 - pd)^2 * exp(-u_d2 * pc * eta);
    S_u0   = 1 - (1 - pd)^2;

    E_us   = 0.5 - 1 / (2 * S_us)   * (1 - pd) * ...
             (exp(-us   * eta * pc * e_det) - exp(-us   * eta * pc * (1 - e_det)));

    E_u_d1 = 0.5 - 1 / (2 * S_u_d1) * (1 - pd) * ...
             (exp(-u_d1 * eta * pc * e_det) - exp(-u_d1 * eta * pc * (1 - e_det)));

    E_u_d2 = 0.5 - 1 / (2 * S_u_d2) * (1 - pd) * ...
             (exp(-u_d2 * eta * pc * e_det) - exp(-u_d2 * eta * pc * (1 - e_det)));

    E_u0   = 0.5;

    S_z = epsilon^2 * S_us + ...
          2 * (1 - epsilon) * epsilon * (1 - (1 - pd)^2 * exp(-us / 2 * pc * eta)) + ...
          (1 - epsilon)^2 * S_u0;

    E_z = (epsilon^2 * S_us + (1 - epsilon)^2 * S_u0) / S_z;

    denom = P_u_d2_2_U * P_u_d1_1_U - P_u_d1_2_L * P_u_d2_1_L;

    s1 = (P_u_d2_2_U * (S_u_d1 - P_u_d1_0_U * S_u0) - ...
          P_u_d1_2_U * (S_u_d2 - P_u_d2_0_L * S_u0)) / denom;

    e1 = (S_u_d1 * E_u_d1 - S_u0 * E_u0 * P_u_d1_0_L) / ...
         (s1 * P_u_d1_1_L);

    R1 = 2 * epsilon * (1 - epsilon) * P_us_1 * s1 * (1 - h(e1)) - ...
         f * S_z * h(E_z);

    if isfinite(R1) && R1 > 0
        R1_all(idx) = R1;
    end
end

% =========================================================================
% 2) 理想光源曲线
% =========================================================================
for idx = 1:numel(d_values)

    d = d_values(idx);

    % 固定 decoy 强度
    decoy_A1_L = 0.01;
    decoy_A2_L = 0.12;
    decoy_B1_L = 0.01;
    decoy_B2_L = 0.12;

    decoy_A1_U = 0.01;
    decoy_A2_U = 0.12;
    decoy_B1_U = 0.01;
    decoy_B2_U = 0.12;

    % ---------------------------------------------------------------------
    % 泊松项
    % ---------------------------------------------------------------------
    P_us_1 = signal_A * exp(-signal_A) * exp(-signal_B) + ...
             signal_B * exp(-signal_B) * exp(-signal_A);

    P_u_d1_0_L = exp(-decoy_A1_L) * exp(-decoy_B1_L);
    P_u_d2_0_L = exp(-decoy_A2_L) * exp(-decoy_B2_L);

    P_u_d1_0_U = exp(-decoy_A1_U) * exp(-decoy_B1_U);
    P_u_d2_0_U = exp(-decoy_A2_U) * exp(-decoy_B2_U);

    P_u_d1_1_L = decoy_A1_L * exp(-decoy_A1_L) * exp(-decoy_B1_L) + ...
                 decoy_B1_L * exp(-decoy_B1_L) * exp(-decoy_A1_L);

    P_u_d2_1_L = decoy_A2_L * exp(-decoy_A2_L) * exp(-decoy_B2_L) + ...
                 decoy_B2_L * exp(-decoy_B2_L) * exp(-decoy_A2_L);

    P_u_d1_1_U = decoy_A1_U * exp(-decoy_A1_U) * exp(-decoy_B1_U) + ...
                 decoy_B1_U * exp(-decoy_B1_U) * exp(-decoy_A1_U);

    P_u_d2_1_U = decoy_A2_U * exp(-decoy_A2_U) * exp(-decoy_B2_U) + ...
                 decoy_B2_U * exp(-decoy_B2_U) * exp(-decoy_A2_U);

    P_u_d1_2_L = exp(-decoy_A1_L) * 0.5 * decoy_B1_L^2 * exp(-decoy_B1_L) + ...
                 exp(-decoy_B1_L) * 0.5 * decoy_A1_L^2 * exp(-decoy_A1_L) + ...
                 decoy_A1_L * exp(-decoy_A1_L) * decoy_B1_L * exp(-decoy_B1_L);

    P_u_d2_2_L = exp(-decoy_A2_L) * 0.5 * decoy_B2_L^2 * exp(-decoy_B2_L) + ...
                 exp(-decoy_B2_L) * 0.5 * decoy_A2_L^2 * exp(-decoy_A2_L) + ...
                 decoy_A2_L * exp(-decoy_A2_L) * decoy_B2_L * exp(-decoy_B2_L);

    P_u_d1_2_U = exp(-decoy_A1_U) * 0.5 * decoy_B1_U^2 * exp(-decoy_B1_U) + ...
                 exp(-decoy_B1_U) * 0.5 * decoy_A1_U^2 * exp(-decoy_A1_U) + ...
                 decoy_A1_U * exp(-decoy_A1_U) * decoy_B1_U * exp(-decoy_B1_U);

    P_u_d2_2_U = exp(-decoy_A2_U) * 0.5 * decoy_B2_U^2 * exp(-decoy_B2_U) + ...
                 exp(-decoy_B2_U) * 0.5 * decoy_A2_U^2 * exp(-decoy_A2_U) + ...
                 decoy_A2_U * exp(-decoy_A2_U) * decoy_B2_U * exp(-decoy_B2_U);

    % ---------------------------------------------------------------------
    % 增益和误码率
    % ---------------------------------------------------------------------
    eta = 10^(-alpha * d / 10);

    us   = signal_A + signal_B;
    u_d1 = decoy_A1_L + decoy_B1_L;
    u_d2 = decoy_A2_L + decoy_B2_L;

    S_us   = 1 - (1 - pd)^2 * exp(-us   * pc * eta);
    S_u_d1 = 1 - (1 - pd)^2 * exp(-u_d1 * pc * eta);
    S_u_d2 = 1 - (1 - pd)^2 * exp(-u_d2 * pc * eta);
    S_u0   = 1 - (1 - pd)^2;

    E_us   = 0.5 - 1 / (2 * S_us)   * (1 - pd) * ...
             (exp(-us   * eta * pc * e_det) - exp(-us   * eta * pc * (1 - e_det)));

    E_u_d1 = 0.5 - 1 / (2 * S_u_d1) * (1 - pd) * ...
             (exp(-u_d1 * eta * pc * e_det) - exp(-u_d1 * eta * pc * (1 - e_det)));

    E_u_d2 = 0.5 - 1 / (2 * S_u_d2) * (1 - pd) * ...
             (exp(-u_d2 * eta * pc * e_det) - exp(-u_d2 * eta * pc * (1 - e_det)));

    E_u0   = 0.5;

    S_z = epsilon^2 * S_us + ...
          2 * (1 - epsilon) * epsilon * (1 - (1 - pd)^2 * exp(-us / 2 * pc * eta)) + ...
          (1 - epsilon)^2 * S_u0;

    E_z = (epsilon^2 * S_us + (1 - epsilon)^2 * S_u0) / S_z;

    denom = P_u_d2_2_U * P_u_d1_1_U - P_u_d1_2_L * P_u_d2_1_L;

    s1 = (P_u_d2_2_U * (S_u_d1 - P_u_d1_0_U * S_u0) - ...
          P_u_d1_2_U * (S_u_d2 - P_u_d2_0_L * S_u0)) / denom;

    e1 = (S_u_d1 * E_u_d1 - S_u0 * E_u0 * P_u_d1_0_L) / ...
         (s1 * P_u_d1_1_L);

    R_ideal = 2 * epsilon * (1 - epsilon) * P_us_1 * s1 * (1 - h(e1)) - ...
              f * S_z * h(E_z);

    if isfinite(R_ideal) && R_ideal > 0
        R_ideal_all(idx) = R_ideal;
    end
end

% =========================================================================
% 3) 同一个 figure 中作图（X 轴放大 2 倍）
% =========================================================================
valid1 = isfinite(R1_all) & (R1_all > 0);
valid2 = isfinite(R_ideal_all) & (R_ideal_all > 0);

x_plot = 2 * d_values;   % 例如 50.5 -> 101

figure;
semilogy(x_plot(valid1), R1_all(valid1), 'r', 'LineWidth', 1.5);
hold on;
semilogy(x_plot(valid2), R_ideal_all(valid2), 'k', 'LineWidth', 1.5);

xlabel('Transmission distance L (km)');
ylabel('Secret key rate R (per pulse)');
title('SNS key rate comparison');
legend('Untrusted source', 'Ideal source', 'Location', 'southwest');
grid on;
xlim([0 1000]);      % 原来 0~500，放大 2 倍后变成 0~1000
ylim([1e-12 1]);