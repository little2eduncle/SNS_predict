function [log_post, d_post_d_thetaR] = logPosterior(thetaP, thetaR, thetaF, C, varargin)
    %logPosterior: Computes the logarithm of the posterior distribution, 
    %              which is the sum of the log-likelihood and the log-prior, 
    %              given the system parameters and observed counts. Optionally 
    %              calculates the gradient of the log-posterior with respect 
    %              to parameters in thetaR.
    %
    % Inputs:
    %   thetaR  - Values of parameters listed in varR 
    %   thetaF  - Values of parameters listed in varF 
    %   varR    - Cell array of random variable names
    %   varF    - Cell array of fixed variable names
    %   C       - Observed click counts
    %   thetaP  - Prior parameters of random variables [alphas, betas, ub, lb]
    %
    % Outputs:
    %   log_posterior - The log-posterior (sum of log-likelihood and log-prior)
    %   d_thetaR      - Derivative of the log-posterior with respect to 
    %                   parameters in thetaR
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).
    % --------------------------------------------------------------------------------------------------------
    % Parse input arguments
    % --------------------------
    p = inputParser;
    addParameter(p, 'fluct_sigma', 0.01);
    parse(p, varargin{:});
    opt = p.Results;

    fluct_sigma      = opt.fluct_sigma;

    if nargout > 1
        % Compute log-likelihood and its derivatives
        [log_likelihood, d_likelihood] = logLikelihood_ZZ(C, thetaR, thetaF, fluct_sigma);
        
        % Compute log-prior and its derivatives
        [log_prior, d_prior] = logPrior(thetaR, thetaP);
        
        % Compute the derivative of the log-posterior

        d_post_d_thetaR = d_likelihood + d_prior;

    else
        % Compute log-likelihood without derivatives
        log_likelihood = logLikelihood_ZZ(C, thetaR, thetaF, fluct_sigma);

        % Compute log-prior without derivatives
        log_prior = logPrior(thetaR, thetaP);
    end

    % Compute the log-posterior as the sum of log-likelihood and log-prior
    log_post = log_likelihood + log_prior;
end