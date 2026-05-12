N_m = 1e7;
frac = 0.043;
pc = 0.6;
pd = 1e-8;
alpha = 0.21;
sig = 0.01;
Nmin=(1-frac)*N_m;
Nmax=(1+frac)*N_m;% 计算[Nmin,Nmax]

delta = UntaggedBitsProb( N_m,sig*N_m,frac );

signal_A = 0.35;
signal_B = 0.35;

decoy_A1 = 0.01;
decoy_A2 = 0.12;

decoy_B1 = 0.01;
decoy_B2 = 0.12;
%----------------------------------------------------------

eta_sa = 2*signal_A/(Nmin+Nmax);
eta_sb = 2*signal_B/(Nmin+Nmax);

eta_d1_a = 2*decoy_A1/(Nmin+Nmax);
eta_d2_a = 2*decoy_A2/(Nmin+Nmax);

eta_d1_b = 2*decoy_B1/(Nmin+Nmax);
eta_d2_b = 2*decoy_B2/(Nmin+Nmax);
%----------------------------------------------------------

N_bar_sa = 1/eta_sa - 1;
N_bar_sb = 1/eta_sb - 1;

N_bar_d1_a = 1/eta_d1_a - 1;
N_bar_d2_a = 1/eta_d2_a - 1;

N_bar_d1_b = 1/eta_d1_b - 1;
N_bar_d2_b = 1/eta_d2_b - 1;
%----------------------------------------------------------

N_bar2_sa = 2/eta_sa - 1;
N_bar2_sb = 2/eta_sb - 1;

N_bar2_d1_a = 2/eta_d1_a - 1;
N_bar2_d2_a = 2/eta_d2_a - 1;

N_bar2_d1_b = 2/eta_d1_b - 1;
N_bar2_d2_b = 2/eta_d2_b - 1;
%----------------------------------------------------------

P_sa_0_L = (1-delta) * (1-eta_sa)^(Nmax);
P_sa_0_U = delta + (1-delta) * (1-eta_sa)^(Nmin) + ...
           delta * (1-eta_sa)^(Nmax+1);

P_sa_1_L = (1-delta) * Nmin * eta_sa * (1-eta_sa)^(Nmin-1);
P_sa_1_U = delta * (Nmin-1) * eta_sa * (1-eta_sa)^(Nmin-2) + ...
           (1-delta) * Nmax * eta_sa * (1-eta_sa)^(Nmax-1) + ...
           delta * N_bar_sa * eta_sa * (1-eta_sa)^(N_bar_sa-1);

P_sa_2_L = (1-delta) * Nmin * (Nmin-1)/2 * eta_sa^2 * (1-eta_sa)^(Nmin-2);
P_sa_2_U = delta * (Nmin-1) * (Nmin-2)/2 * eta_sa^2 * (1-eta_sa)^(Nmin-3) +...
           (1-delta) * Nmax * (Nmax-1)/2 * eta_sa^2 * (1-eta_sa)^(Nmax-2) +...
           delta * N_bar2_sa * (N_bar2_sa-1)/2 * eta_sa^2 * (1-eta_sa)^(N_bar2_sa-2);

sa_L = [P_sa_0_L, P_sa_1_L, P_sa_2_L];
sa_U = [P_sa_0_U, P_sa_1_U, P_sa_2_U];
%----------------------------------------------------------

P_d1_a_0_L = (1-delta) * (1-eta_d1_a)^(Nmax);
P_d1_a_0_U = delta + (1-delta) * (1-eta_d1_a)^(Nmin) + ...
           delta * (1-eta_d1_a)^(Nmax+1);

P_d1_a_1_L = (1-delta) * Nmin * eta_d1_a * (1-eta_d1_a)^(Nmin-1);
P_d1_a_1_U = delta * (Nmin-1) * eta_d1_a * (1-eta_d1_a)^(Nmin-2) + ...
           (1-delta) * Nmax * eta_d1_a * (1-eta_d1_a)^(Nmax-1) + ...
           delta * N_bar_d1_a * eta_d1_a * (1-eta_d1_a)^(N_bar_d1_a-1);

P_d1_a_2_L = (1-delta) * Nmin * (Nmin-1)/2 * eta_d1_a^2 * (1-eta_d1_a)^(Nmin-2);
P_d1_a_2_U = delta * (Nmin-1) * (Nmin-2)/2 * eta_d1_a^2 * (1-eta_d1_a)^(Nmin-3) +...
           (1-delta) * Nmax * (Nmax-1)/2 * eta_d1_a^2 * (1-eta_d1_a)^(Nmax-2) +...
           delta * N_bar2_d1_a * (N_bar2_d1_a-1)/2 * eta_d1_a^2 * (1-eta_d1_a)^(N_bar2_d1_a-2);

d1_a_L = [P_d1_a_0_L, P_d1_a_1_L, P_d1_a_2_L];
d1_a_U = [P_d1_a_0_U, P_d1_a_1_U, P_d1_a_2_U];
%----------------------------------------------------------

P_d2_a_0_L = (1-delta) * (1-eta_d2_a)^(Nmax);
P_d2_a_0_U = delta + (1-delta) * (1-eta_d2_a)^(Nmin) + ...
           delta * (1-eta_d2_a)^(Nmax+1);

P_d2_a_1_L = (1-delta) * Nmin * eta_d2_a * (1-eta_d2_a)^(Nmin-1);
P_d2_a_1_U = delta * (Nmin-1) * eta_d2_a * (1-eta_d2_a)^(Nmin-2) + ...
           (1-delta) * Nmax * eta_d2_a * (1-eta_d2_a)^(Nmax-1) + ...
           delta * N_bar_d2_a * eta_d2_a * (1-eta_d2_a)^(N_bar_d2_a-1);

P_d2_a_2_L = (1-delta) * Nmin * (Nmin-1)/2 * eta_d2_a^2 * (1-eta_d2_a)^(Nmin-2);
P_d2_a_2_U = delta * (Nmin-1) * (Nmin-2)/2 * eta_d2_a^2 * (1-eta_d2_a)^(Nmin-3) +...
           (1-delta) * Nmax * (Nmax-1)/2 * eta_d2_a^2 * (1-eta_d2_a)^(Nmax-2) +...
           delta * N_bar2_d2_a * (N_bar2_d2_a-1)/2 * eta_d2_a^2 * (1-eta_d2_a)^(N_bar2_d2_a-2);

d2_a_L = [P_d2_a_0_L, P_d2_a_1_L, P_d2_a_2_L];
d2_a_U = [P_d2_a_0_U, P_d2_a_1_U, P_d2_a_2_U];
%----------------------------------------------------------
%----------------------------------------------------------

P_sb_0_L = (1-delta) * (1-eta_sb)^(Nmax);
P_sb_0_U = delta + (1-delta) * (1-eta_sb)^(Nmin) + ...
           delta * (1-eta_sb)^(Nmax+1);

P_sb_1_L = (1-delta) * Nmin * eta_sb * (1-eta_sb)^(Nmin-1);
P_sb_1_U = delta * (Nmin-1) * eta_sb * (1-eta_sb)^(Nmin-2) + ...
           (1-delta) * Nmax * eta_sb * (1-eta_sb)^(Nmax-1) + ...
           delta * N_bar_sb * eta_sb * (1-eta_sb)^(N_bar_sb-1);

P_sb_2_L = (1-delta) * Nmin * (Nmin-1)/2 * eta_sb^2 * (1-eta_sb)^(Nmin-2);
P_sb_2_U = delta * (Nmin-1) * (Nmin-2)/2 * eta_sb^2 * (1-eta_sb)^(Nmin-3) +...
           (1-delta) * Nmax * (Nmax-1)/2 * eta_sb^2 * (1-eta_sb)^(Nmax-2) +...
           delta * N_bar2_sb * (N_bar2_sb-1)/2 * eta_sb^2 * (1-eta_sb)^(N_bar2_sb-2);

sb_L = [P_sb_0_L, P_sb_1_L, P_sb_2_L];
sb_U = [P_sb_0_U, P_sb_1_U, P_sb_2_U];
%----------------------------------------------------------

P_d1_b_0_L = (1-delta) * (1-eta_d1_b)^(Nmax);
P_d1_b_0_U = delta + (1-delta) * (1-eta_d1_b)^(Nmin) + ...
           delta * (1-eta_d1_b)^(Nmax+1);

P_d1_b_1_L = (1-delta) * Nmin * eta_d1_b * (1-eta_d1_b)^(Nmin-1);
P_d1_b_1_U = delta * (Nmin-1) * eta_d1_b * (1-eta_d1_b)^(Nmin-2) + ...
           (1-delta) * Nmax * eta_d1_b * (1-eta_d1_b)^(Nmax-1) + ...
           delta * N_bar_d1_b * eta_d1_b * (1-eta_d1_b)^(N_bar_d1_b-1);

P_d1_b_2_L = (1-delta) * Nmin * (Nmin-1)/2 * eta_d1_b^2 * (1-eta_d1_b)^(Nmin-2);
P_d1_b_2_U = delta * (Nmin-1) * (Nmin-2)/2 * eta_d1_b^2 * (1-eta_d1_b)^(Nmin-3) +...
           (1-delta) * Nmax * (Nmax-1)/2 * eta_d1_b^2 * (1-eta_d1_b)^(Nmax-2) +...
           delta * N_bar2_d1_b * (N_bar2_d1_b-1)/2 * eta_d1_b^2 * (1-eta_d1_b)^(N_bar2_d1_b-2);

d1_b_L = [P_d1_b_0_L, P_d1_b_1_L, P_d1_b_2_L];
d1_b_U = [P_d1_b_0_U, P_d1_b_1_U, P_d1_b_2_U];
%----------------------------------------------------------

P_d2_b_0_L = (1-delta) * (1-eta_d2_b)^(Nmax);
P_d2_b_0_U = delta + (1-delta) * (1-eta_d2_b)^(Nmin) + ...
           delta * (1-eta_d2_b)^(Nmax+1);

P_d2_b_1_L = (1-delta) * Nmin * eta_d2_b * (1-eta_d2_b)^(Nmin-1);
P_d2_b_1_U = delta * (Nmin-1) * eta_d2_b * (1-eta_d2_b)^(Nmin-2) + ...
           (1-delta) * Nmax * eta_d2_b * (1-eta_d2_b)^(Nmax-1) + ...
           delta * N_bar_d2_b * eta_d2_b * (1-eta_d2_b)^(N_bar_d2_b-1);

P_d2_b_2_L = (1-delta) * Nmin * (Nmin-1)/2 * eta_d2_b^2 * (1-eta_d2_b)^(Nmin-2);
P_d2_b_2_U = delta * (Nmin-1) * (Nmin-2)/2 * eta_d2_b^2 * (1-eta_d2_b)^(Nmin-3) +...
           (1-delta) * Nmax * (Nmax-1)/2 * eta_d2_b^2 * (1-eta_d2_b)^(Nmax-2) +...
           delta * N_bar2_d2_b * (N_bar2_d2_b-1)/2 * eta_d2_b^2 * (1-eta_d2_b)^(N_bar2_d2_b-2);

d2_b_L = [P_d2_b_0_L, P_d2_b_1_L, P_d2_b_2_L];
d2_b_U = [P_d2_b_0_U, P_d2_b_1_U, P_d2_b_2_U];
%----------------------------------------------------------

S = load('thetaR_ZZ_MAP_results.mat');

x = zeros(1, length(S.L_values));

for i = 1:length(S.L_values)

    L = S.L_values(i);
    theta = S.thetaR_MAP_mat(i,:);

    if(theta(1)<=0)
        break;
    end
    signal_A = theta(1);
    signal_B = theta(2);

    x(i) = signal_A*exp(-signal_A)*exp(-signal_B) + ...
           signal_B*exp(-signal_B)*exp(-signal_A);

    % fprintf('L = %.1f km: A = %.4f, B = %.4f\n', L, signal_A, signal_B);

end

P_s_1_L = P_sa_1_L*P_sb_0_L + P_sa_0_L*P_sb_1_L;
P_s_1_U = P_sa_1_U*P_sb_0_U + P_sa_0_U*P_sb_1_U;

Low = ones(1, length(S.L_values));
Low = P_s_1_L * Low;

Up = ones(1, length(S.L_values));
Up = P_s_1_U * Up;

step = 5;

idx_valid = (x > 1e-12);

L_plot = S.L_values(idx_valid);
x_plot = x(idx_valid);
Low_plot = Low(idx_valid);
Up_plot = Up(idx_valid);

% 下采样
L_plot = L_plot(1:step:end);
x_plot = x_plot(1:step:end);
Low_plot = Low_plot(1:step:end);
Up_plot = Up_plot(1:step:end);

figure;
plot(L_plot, x_plot, 'o-', 'LineWidth', 1.5); hold on;
plot(L_plot, Low_plot, '--', 'LineWidth', 1.5);
plot(L_plot, Up_plot, '--', 'LineWidth', 1.5);

xlabel('Distance (km)');
ylabel('Value');
legend('x (MAP)', 'Lower Bound', 'Upper Bound');
grid on;