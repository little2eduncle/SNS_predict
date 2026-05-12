function [log_gamma, d_gamma] = logGamma(x, alpha, beta)
    %logGamma: Computes the logarithm of the Gamma distribution for a given 
    %          input x and parameters alpha (shape) and beta (rate). Optionally 
    %          calculates the derivative of the log-Gamma w.r.t the parameters.
    %
    % Inputs:
    %   x      - Input value(s) where the Gamma distribution is evaluated
    %   alpha  - Shape parameter of the Gamma distribution
    %   beta   - Rate parameter of the Gamma distribution
    %
    % Outputs:
    %   log_gamm - Logarithm of the Gamma distribution evaluated at x

    % Adjust x to ensure numerical stability (ensures x > 0)
    x = x + eps;

    % Compute the logarithm of the Gamma distribution (without normalization)
    log_gamma = (alpha - 1) .* log(x + (alpha == 1)) - beta .* x; % + alpha * log(beta) - gammaln(alpha)
    % Derivative of log-Gamma w.r.t. x
    d_gamma = (alpha - 1) ./ (x + (alpha == 1)) - beta;
end
