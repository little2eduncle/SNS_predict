function [log_likelihood_ZZ, d_log_likelihood_ZZ] = logLikelihood_ZZ(Nmin, Nmax, thetaR_ZZ)

    [signal_A, signal_B] = deal(thetaR_ZZ{:});

    S = Nmin + Nmax;

    eta_signal_A = 2 * signal_A / S;
    eta_signal_B = 2 * signal_B / S;

    log_likelihood_ZZ = Nmax .* log(1 - eta_signal_A) + log(eta_signal_B) + (Nmin - 1) .* log(1 - eta_signal_B) + ...
                        Nmax .* log(1 - eta_signal_B) + log(eta_signal_A) + (Nmin - 1) .* log(1 - eta_signal_A);

    d_log_likelihood_ZZ = zeros(1, 2);

    d_log_likelihood_ZZ(1) = (2 / S) * ( ...
        1 ./ eta_signal_A - (S - 1) ./ (1 - eta_signal_A) );

    d_log_likelihood_ZZ(2) = (2 / S) * ( ...
        1 ./ eta_signal_B - (S - 1) ./ (1 - eta_signal_B) );
end