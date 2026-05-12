function compare_simulation_with_fluct_analytic()
% Compare simulate_SNS_TF_QKD_3 against
%   Zwindow_click_prob_asym_fluct
%   XXwindow_click_prob_asym_screened_fluct
%
% Make sure these functions are on MATLAB path:
%   simulate_SNS_TF_QKD_3
%   Zwindow_click_prob_asym_fluct
%   XXwindow_click_prob_asym_screened_fluct

    clc;
    rng(1);

    % -------------------------------------------------------------
    % 1) User parameters
    % -------------------------------------------------------------
    N = 1e7;              % increase if you want smaller Monte Carlo noise
    nGH = 15;             % Gauss-Hermite order for fluctuation averaging

    signal_A = 0.35;
    signal_B = 0.35;

    decoy_A  = [0.01, 0.12];
    decoy_B  = [0.01, 0.12];

    alpha_A  = 0.2;
    alpha_B  = 0.2;
    dAC      = 50;
    dBC      = 50;

    pa0 = 0.01;
    pa1 = 0.01;
    pc0 = 0.20;
    pc1 = 0.22;
    pd0 = 1e-6;
    pd1 = 1e-6;

    z_window    = 0.5;
    epsilon     = 0.99;
    lambda      = 0.2;
    fluct_sigma = 0.02;
    vacuum      = 0;

    % X-window probabilities for [decoy_1, decoy_2, ..., vacuum]
    p_decoy = [0.45, 0.45, 0.10];

    thetaA = {signal_A, decoy_A, alpha_A, dAC};
    thetaB = {signal_B, decoy_B, alpha_B, dBC};
    thetaC = {pa0, pa1, pc0, pc1, pd0, pd1};

    % -------------------------------------------------------------
    % 2) Run simulation
    % -------------------------------------------------------------
    [C_total, C_windows, C_ZZ_intensity, C_XX_intensity] = ...
        simulate_SNS_TF_QKD_4( ...
            N, thetaA, thetaB, thetaC, ...
            'z_window', z_window, ...
            'epsilon', epsilon, ...
            'p_decoy', p_decoy, ...
            'lambda', lambda, ...
            'fluct_sigma', fluct_sigma, ...
            'vacuum', vacuum, ...
            'use_gpu', false);

    fprintf('Simulation finished.\n');
    fprintf('C_total = [%g %g %g %g]\n', C_total);
    fprintf('C_windows.ZZ = [%g %g %g %g]\n', C_windows.ZZ);
    fprintf('C_windows.XX = [%g %g %g %g]\n\n', C_windows.XX);

    % -------------------------------------------------------------
    % 3) Compare ZZ window
    % -------------------------------------------------------------
    fprintf('==================== ZZ comparison ====================\n');

    muZZ_A = C_ZZ_intensity.mu_nominal_A(:);
    muZZ_B = C_ZZ_intensity.mu_nominal_B(:).';

    for i = 1:numel(muZZ_A)
        for j = 1:numel(muZZ_B)
            muA0 = muZZ_A(i);
            muB0 = muZZ_B(j);

            total_valid = C_ZZ_intensity.total_valid(i,j);

            sim_D1 = safe_div(C_ZZ_intensity.D1_only(i,j), total_valid);
            sim_D0 = safe_div(C_ZZ_intensity.D0_only(i,j), total_valid);
            sim_rest = safe_div(C_ZZ_intensity.rest(i,j), total_valid);

            [rZ, ~] = Zwindow_click_prob_asym_fluct( ...
                muA0, muB0, ...
                pd0, pc0, pa0, ...
                pd1, pc1, pa1, ...
                alpha_A, alpha_B, dAC, dBC, ...
                fluct_sigma, ...
                'nGH', nGH);

            fprintf('ZZ pair (muA=%.4f, muB=%.4f), total=%g\n', muA0, muB0, total_valid);
            fprintf('  sim  : D0=%.8f, D1=%.8f, rest=%.8f\n', sim_D0, sim_D1, sim_rest);
            fprintf('  anal : D0=%.8f, D1=%.8f, rest=%.8f\n', ...
                rZ.P_only_D0, rZ.P_only_D1, rZ.P_rest);
            fprintf('  diff : D0=%+.3e, D1=%+.3e, rest=%+.3e\n\n', ...
                sim_D0 - rZ.P_only_D0, ...
                sim_D1 - rZ.P_only_D1, ...
                sim_rest - rZ.P_rest);
        end
    end

    % -------------------------------------------------------------
    % 4) Compare XX window
    % -------------------------------------------------------------
    fprintf('==================== XX comparison ====================\n');

    muXX_A = C_XX_intensity.mu_nominal_A(:);
    muXX_B = C_XX_intensity.mu_nominal_B(:).';

    for i = 1:numel(muXX_A)
        for j = 1:numel(muXX_B)
            muA0 = muXX_A(i);
            muB0 = muXX_B(j);

            total_valid = C_XX_intensity.total_valid(i,j);
            total_raw   = C_XX_intensity.total_raw(i,j);

            % These are CONDITIONAL on accepted phase
            sim_D1 = safe_div(C_XX_intensity.D1_only(i,j), total_valid);
            sim_D0 = safe_div(C_XX_intensity.D0_only(i,j), total_valid);
            sim_rest = safe_div(C_XX_intensity.rest(i,j), total_valid);

            % This is acceptance fraction for this XX pair
            sim_accept = safe_div(total_valid, total_raw);

            [rXX, ~] = XXwindow_click_prob_asym_screened_fluct( ...
                muA0, muB0, ...
                pd0, pc0, pa0, ...
                pd1, pc1, pa1, ...
                alpha_A, alpha_B, dAC, dBC, ...
                lambda, fluct_sigma, ...
                'nGH', nGH);

            fprintf('XX pair (muA=%.4f, muB=%.4f), raw=%g, valid=%g\n', ...
                muA0, muB0, total_raw, total_valid);

            fprintf('  sim(cond)  : D0=%.8f, D1=%.8f, rest=%.8f\n', sim_D0, sim_D1, sim_rest);
            fprintf('  anal(cond) : D0=%.8f, D1=%.8f, rest=%.8f\n', ...
                rXX.P_only_D0, rXX.P_only_D1, rXX.P_rest);
            fprintf('  diff(cond) : D0=%+.3e, D1=%+.3e, rest=%+.3e\n', ...
                sim_D0 - rXX.P_only_D0, ...
                sim_D1 - rXX.P_only_D1, ...
                sim_rest - rXX.P_rest);

            fprintf('  sim accept = %.8f\n', sim_accept);
            fprintf('  anal accept= %.8f\n', rXX.acceptFraction);
            fprintf('  diff accept= %+.3e\n\n', sim_accept - rXX.acceptFraction);
        end
    end
end


function y = safe_div(a, b)
    if b <= 0
        y = NaN;
    else
        y = a / b;
    end
end