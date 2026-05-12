function [s1, e1, P] = Z1_data_UB(signal_A, decoy_A, signal_B, decoy_B, N_m, frac, sig, pd, pc, alpha, d, e_det)

Nmin=(1-frac)*N_m;
Nmax=(1+frac)*N_m;% 计算[Nmin,Nmax]

eta_trans = 10^(-alpha*d/10)*pc;

delta = UntaggedBitsProb( N_m,sig*N_m,frac );

%----------------------------------------------------------

eta_sa = 2*signal_A/(Nmin+Nmax);
eta_sb = 2*signal_B/(Nmin+Nmax);

eta_d1_a = 2*decoy_A(1)/(Nmin+Nmax);
eta_d2_a = 2*decoy_A(2)/(Nmin+Nmax);

eta_d1_b = 2*decoy_B(1)/(Nmin+Nmax);
eta_d2_b = 2*decoy_B(2)/(Nmin+Nmax);
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

%----------------------------------------------------------

P_s_0_L = P_sa_0_L * P_sb_0_L;
P_s_0_U = P_sa_0_U * P_sb_0_U;

P_s_1_L = P_sa_0_L * P_sb_1_L + P_sa_1_L * P_sb_0_L;
P_s_1_U = P_sa_0_U * P_sb_1_U + P_sa_1_U * P_sb_0_U;

P_s_2_L = P_sa_0_L * P_sb_2_L + P_sa_2_L * P_sb_0_L + ...
          P_sa_1_L * P_sb_1_L;
P_s_2_U = P_sa_0_U * P_sb_2_U + P_sa_2_U * P_sb_0_U + ...
          P_sa_1_U * P_sb_1_U;

P_s_1 = [P_s_1_L, P_s_1_U];
%----------------------------------------------------------

P_d1_0_L = P_d1_a_0_L * P_d1_b_0_L;
P_d1_0_U = P_d1_a_0_U * P_d1_b_0_U;

P_d1_1_L = P_d1_a_0_L * P_d1_b_1_L + P_d1_a_1_L * P_d1_b_0_L;
P_d1_1_U = P_d1_a_0_U * P_d1_b_1_U + P_d1_a_1_U * P_d1_b_0_U;

P_d1_2_L = P_d1_a_0_L * P_d1_b_2_L + P_d1_a_2_L * P_d1_b_0_L + ...
          P_d1_a_1_L * P_d1_b_1_L;
P_d1_2_U = P_d1_a_0_U * P_d1_b_2_U + P_d1_a_2_U * P_d1_b_0_U + ...
          P_d1_a_1_U * P_d1_b_1_U;

P_d1_1 = [P_d1_1_L, P_d1_1_U];
P_d1_0 = [P_d1_0_L, P_d1_0_U];
P_d1_2 = [P_d1_2_L, P_d1_2_U];
%----------------------------------------------------------

P_d2_0_L = P_d2_a_0_L * P_d2_b_0_L;
P_d2_0_U = P_d2_a_0_U * P_d2_b_0_U;

P_d2_1_L = P_d2_a_0_L * P_d2_b_1_L + P_d2_a_1_L * P_d2_b_0_L;
P_d2_1_U = P_d2_a_0_U * P_d2_b_1_U + P_d2_a_1_U * P_d2_b_0_U;

P_d2_2_L = P_d2_a_0_L * P_d2_b_2_L + P_d2_a_2_L * P_d2_b_0_L + ...
          P_d2_a_1_L * P_d2_b_1_L;
P_d2_2_U = P_d2_a_0_U * P_d2_b_2_U + P_d2_a_2_U * P_d2_b_0_U + ...
          P_d2_a_1_U * P_d2_b_1_U;

P_d2_0 = [P_d2_0_L, P_d2_0_U];
P_d2_1 = [P_d2_1_L, P_d2_1_U];
P_d2_2 = [P_d2_2_L, P_d2_2_U];
%----------------------------------------------------------
S_s = 1 - (1-pd)^2*exp(-2*signal_A*eta_trans);
S_s_half = 1 - (1-pd)^2*exp(-signal_A*eta_trans);
S_d1 = 1 - (1-pd)^2*exp(-2*decoy_A(1)*eta_trans);
S_d2 = 1 - (1-pd)^2*exp(-2*decoy_A(2)*eta_trans);
S_0 = 1 - (1-pd)^2;

%----------------------------------------------------------
E_s = 0.5 - 1/(2*S_s) * (1-pd) * ...
     (exp(-2*signal_A*eta_trans*e_det) - ...
     exp(-2*signal_A*eta_trans*(1-e_det)));

E_d1 = 0.5 - 1/(2*S_d1) * (1-pd) * ...
     (exp(-2*decoy_A(1)*eta_trans*e_det) - ...
     exp(-2*decoy_A(1)*eta_trans*(1-e_det)));

E_d2 = 0.5 - 1/(2*S_d2) * (1-pd) * ...
     (exp(-2*decoy_A(2)*eta_trans*e_det) - ...
     exp(-2*decoy_A(2)*eta_trans*(1-e_det)));

E_0 = 0.5;

%----------------------------------------------------------
s1 = (P_d2_2_L * (S_d1 - P_d1_0_U * S_0) - P_d1_2_U * (S_d2 - P_d2_0_L * S_0)) / ...
      (P_d2_2_U * P_d1_1_U - P_d1_2_L * P_d2_1_L);

e1 = (S_d1 * E_d1 - S_0 * E_0 * P_d1_0_L) /...
      (s1 * P_d1_1_L);

P = {P_s_1, P_d1_1, P_d1_0, P_d1_2, P_d2_0, P_d2_1, P_d2_2};
%----------------------------------------------------------
end