function [s1, e1, P] = Z1_data_Bayes(signal_A, decoy_A, signal_B, decoy_B, pd, pc, alpha, d, e_det)
sA = signal_A;
sB = signal_B;

A1 = decoy_A(1);
A2 = decoy_A(2);

B1 = decoy_B(1);
B2 = decoy_B(2);

signal_expect = 0.496;
decoy1_expect = 0.1;
decoy2_expect = 0.411;

eta_trans = 10^(-alpha*d/10)*pc;
%----------------------------------------------------------

P_sa_0_L = exp(-sA);
P_sa_0_U = exp(-sA);

P_sa_1_L = sA*exp(-sA);
P_sa_1_U = sA*exp(-sA);

P_sa_2_L = 0.5*sA^2*exp(-sA);
P_sa_2_U = 0.5*sA^2*exp(-sA);
%----------------------------------------------------------

P_d1_a_0_L = exp(-A1);
P_d1_a_0_U = exp(-A1);

P_d1_a_1_L = A1*exp(-A1);
P_d1_a_1_U = A1*exp(-A1);

P_d1_a_2_L = 0.5*A1^2*exp(-A1);
P_d1_a_2_U = 0.5*A1^2*exp(-A1);
%----------------------------------------------------------

P_d2_a_0_L = exp(-A2);
P_d2_a_0_U = exp(-A2);

P_d2_a_1_L = A2*exp(-A2);
P_d2_a_1_U = A2*exp(-A2);

P_d2_a_2_L = 0.5*A2^2*exp(-A2);
P_d2_a_2_U = 0.5*A2^2*exp(-A2);
%----------------------------------------------------------

P_sb_0_L = exp(-sB);
P_sb_0_U = exp(-sB);

P_sb_1_L = sB*exp(-sB);
P_sb_1_U = sB*exp(-sB);

P_sb_2_L = 0.5*sB^2*exp(-sB);
P_sb_2_U = 0.5*sB^2*exp(-sB);
%----------------------------------------------------------

P_d1_b_0_L = exp(-B1);
P_d1_b_0_U = exp(-B1);

P_d1_b_1_L = B1*exp(-B1);
P_d1_b_1_U = B1*exp(-B1);

P_d1_b_2_L = 0.5*B1^2*exp(-B1);
P_d1_b_2_U = 0.5*B1^2*exp(-B1);
%----------------------------------------------------------

P_d2_b_0_L = exp(-B2);
P_d2_b_0_U = exp(-B2);

P_d2_b_1_L = B2*exp(-B2);
P_d2_b_1_U = B2*exp(-B2);

P_d2_b_2_L = 0.5*B2^2*exp(-B2);
P_d2_b_2_U = 0.5*B2^2*exp(-B2);

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
S_s = 1 - (1-pd)^2*exp(-2*signal_expect*eta_trans);
S_s_half = 1 - (1-pd)^2*exp(-signal_expect*eta_trans);
S_d1 = 1 - (1-pd)^2*exp(-2*decoy1_expect*eta_trans);
S_d2 = 1 - (1-pd)^2*exp(-2*decoy2_expect*eta_trans);
S_0 = 1 - (1-pd)^2;

%----------------------------------------------------------
E_s = 0.5 - 1/(2*S_s) * (1-pd) * ...
     (exp(-2*signal_expect*eta_trans*e_det) - ...
     exp(-2*signal_expect*eta_trans*(1-e_det)));

E_d1 = 0.5 - 1/(2*S_d1) * (1-pd) * ...
     (exp(-2*decoy1_expect*eta_trans*e_det) - ...
     exp(-2*decoy1_expect*eta_trans*(1-e_det)));

E_d2 = 0.5 - 1/(2*S_d2) * (1-pd) * ...
     (exp(-2*decoy2_expect*eta_trans*e_det) - ...
     exp(-2*decoy2_expect*eta_trans*(1-e_det)));

E_0 = 0.5;

%----------------------------------------------------------
s1 = (P_d2_2_L * (S_d1 - P_d1_0_U * S_0) - P_d1_2_U * (S_d2 - P_d2_0_L * S_0)) / ...
      (P_d2_2_U * P_d1_1_U - P_d1_2_L * P_d2_1_L);

e1 = (S_d1 * E_d1 - S_0 * E_0 * P_d1_0_L) /...
      (s1 * P_d1_1_L);

P = {P_s_1, P_d1_1, P_d1_0, P_d1_2, P_d2_0, P_d2_1, P_d2_2};
end