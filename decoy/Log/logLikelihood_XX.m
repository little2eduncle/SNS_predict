function [log_likelihood_XX, d_log_likelihood_XX] = logLikelihood_XX(C, thetaR_XX, thetaF_XX, lambda, fluct_sigma)

    [decoy_A, decoy_B] = deal(thetaR_XX{:});
    [alpha_A, alpha_B, dAC, dBC, pa0, pa1, pc0, pc1, pd0, pd1, ~, ~] = deal(thetaF_XX{:});

    decoy_A = decoy_A(:).';
    decoy_B = decoy_B(:).';

    nA = numel(decoy_A);
    nB = numel(decoy_B);

    log_likelihood_XX = 0;
    d_log_likelihood_XX = zeros(1, nA + nB);

    epsP = 1e-300;

    for i = 1:nA

        [result, grad] = XXwindow_click_prob_asym_screened_fluct( ...
                decoy_A(i), decoy_B(i), ...
                pd0, pc0, pa0, ...
                pd1, pc1, pa1, ...
                alpha_A, alpha_B, dAC, dBC, ...
                lambda, fluct_sigma, ...
                'nGH', 15);

            % Model probability order:
            % [only D0, only D1, rest]
            P3 = [result.P_only_D0, result.P_only_D1, result.P_rest];
            P3 = max(P3, epsP);

            % Simulation event order:
            % [D1_only, D0_only, rest]
            Ck = squeeze(C(i,i,:)).';

            % Compress to:
            % [only D0, only D1, rest]
            C_counts = [Ck(2), Ck(1), Ck(3)];

            log_likelihood_XX = log_likelihood_XX + sum(C_counts .* log(P3));

            % grad.muA = [dP_only_D0, dP_only_D1, dP_rest]
            % grad.muB = [dP_only_D0, dP_only_D1, dP_rest]
            
            d_log_likelihood_XX(i) = d_log_likelihood_XX(i)    + sum(C_counts .* grad.muA ./ P3);

            d_log_likelihood_XX(nA+i) = d_log_likelihood_XX(nA+i) + sum(C_counts .* grad.muB ./ P3);
    end

end