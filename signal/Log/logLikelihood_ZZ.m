function [log_likelihood_ZZ, d_log_likelihood_ZZ] = logLikelihood_ZZ(C_ZZ, thetaR_ZZ, thetaF_ZZ, fluct_sigma)

    [signal_A, signal_B] = deal(thetaR_ZZ{:});
    [alpha_A, alpha_B, dAC, dBC, ~, ~, pa0, pa1, pc0, pc1, pd0, pd1] = deal(thetaF_ZZ{:});

    d_log_likelihood_ZZ = zeros(1, 2);

    [result_A, grad_A] = Zwindow_click_prob_asym_fluct( ...
        signal_A, 0, ...
        pd0, pc0, pa0, ...
        pd1, pc1, pa1, ...
        alpha_A, alpha_B, dAC, dBC, ...
        fluct_sigma, ...
        'nGH', 15);

    [result_B, grad_B] = Zwindow_click_prob_asym_fluct( ...
        0, signal_B, ...
        pd0, pc0, pa0, ...
        pd1, pc1, pa1, ...
        alpha_A, alpha_B, dAC, dBC, ...
        fluct_sigma, ...
        'nGH', 15);

    % Model probability order:
    % [only D0, only D1, rest]
    P3_A = [result_A.P_only_D0, result_A.P_only_D1, result_A.P_rest];
    P3_B = [result_B.P_only_D0, result_B.P_only_D1, result_B.P_rest];

    P3_A = max(P3_A, eps);
    P3_B = max(P3_B, eps);

    % Simulation event order:
    % [D1_only, D0_only, rest]
    Ck_A = squeeze(C_ZZ(1,2,:)).';
    Ck_B = squeeze(C_ZZ(2,1,:)).';

    % Compress to:
    % [only D0, only D1, rest]
    CA_counts = [Ck_A(2), Ck_A(1), Ck_A(3)];
    CB_counts = [Ck_B(2), Ck_B(1), Ck_B(3)];

    log_likelihood_ZZ = sum(CA_counts .* log(P3_A)) + sum(CB_counts .* log(P3_B));
    
    d_log_likelihood_ZZ(1) = sum(CA_counts .* grad_A.muA ./ P3_A);
    d_log_likelihood_ZZ(2) = sum(CB_counts .* grad_B.muB ./ P3_B);

end