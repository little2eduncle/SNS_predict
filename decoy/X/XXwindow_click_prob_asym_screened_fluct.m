function [result, grad] = XXwindow_click_prob_asym_screened_fluct( ...
    muA0, muB0, ...
    pd0, pc0, pa0, ...
    pd1, pc1, pa1, ...
    alpha_A, alpha_B, dAC, dBC, ...
    lambda, fluct_sigma, varargin)
% XXwindow_click_prob_asym_screened_fluct
%
% Fluctuation-averaged version of XXwindow_click_prob_asym_screened.
%
% The fluctuation model matches the simulator:
%   mu = max(mu0 + N(0, fluct_sigma^2), 0) for nonzero nominal intensities
%   vacuum (mu0 = 0) does not fluctuate
%
% Optional:
%   'nGH' : Gauss-Hermite quadrature order, default 15

    ip = inputParser;
    addParameter(ip, 'nGH', 15);
    parse(ip, varargin{:});
    nGH = ip.Results.nGH;

    assert(muA0 >= 0 && muB0 >= 0, 'muA0 and muB0 must be nonnegative.');
    assert(fluct_sigma >= 0, 'fluct_sigma must be nonnegative.');
    assert(lambda > 0 && lambda <= 1, 'lambda must satisfy 0 < lambda <= 1.');

    [muA_nodes, wA, dmuA_dmu0, dmuA_dsigma] = gaussian_fluct_rule(muA0, fluct_sigma, nGH);
    [muB_nodes, wB, dmuB_dmu0, dmuB_dsigma] = gaussian_fluct_rule(muB0, fluct_sigma, nGH);
    
    result = init_result_XX_fluct();
    grad   = init_grad_rowvec_fluct_XX();

    tauA = 10^(-alpha_A * dAC / 10);
    tauB = 10^(-alpha_B * dBC / 10);

    mean_muA = 0;
    mean_muB = 0;

    for i = 1:numel(muA_nodes)
        for j = 1:numel(muB_nodes)
            wt = wA(i) * wB(j);

            muAi = muA_nodes(i);
            muBj = muB_nodes(j);

            [r, g] = XXwindow_click_prob_asym_screened( ...
                muAi, muBj, ...
                pd0, pc0, pa0, ...
                pd1, pc1, pa1, ...
                alpha_A, alpha_B, dAC, dBC, ...
                lambda);

            % ---- conditional accepted probabilities
            result.P_only_D0 = result.P_only_D0 + wt * r.P_only_D0;
            result.P_only_D1 = result.P_only_D1 + wt * r.P_only_D1;
            result.P_none    = result.P_none    + wt * r.P_none;
            result.P_double  = result.P_double  + wt * r.P_double;
            result.P_rest    = result.P_rest    + wt * r.P_rest;

            % ---- accepted unconditional probabilities
            result.P_only_D0_valid = result.P_only_D0_valid + wt * r.P_only_D0_valid;
            result.P_only_D1_valid = result.P_only_D1_valid + wt * r.P_only_D1_valid;
            result.P_none_valid    = result.P_none_valid    + wt * r.P_none_valid;
            result.P_double_valid  = result.P_double_valid  + wt * r.P_double_valid;
            result.P_rest_valid    = result.P_rest_valid    + wt * r.P_rest_valid;

            % ---- gradients wrt nominal muA/muB
            grad.muA = grad.muA + wt * g.muA * dmuA_dmu0(i);
            grad.muB = grad.muB + wt * g.muB * dmuB_dmu0(j);

            % ---- gradient wrt fluct_sigma
            grad.fluct_sigma = grad.fluct_sigma + wt * ( ...
                g.muA * dmuA_dsigma(i) + g.muB * dmuB_dsigma(j));

            % ---- direct passthrough gradients
            grad.pd0     = grad.pd0     + wt * g.pd0;
            grad.pc0     = grad.pc0     + wt * g.pc0;
            grad.pa0     = grad.pa0     + wt * g.pa0;
            grad.pd1     = grad.pd1     + wt * g.pd1;
            grad.pc1     = grad.pc1     + wt * g.pc1;
            grad.pa1     = grad.pa1     + wt * g.pa1;
            grad.alpha_A = grad.alpha_A + wt * g.alpha_A;
            grad.alpha_B = grad.alpha_B + wt * g.alpha_B;
            grad.dAC     = grad.dAC     + wt * g.dAC;
            grad.dBC     = grad.dBC     + wt * g.dBC;
            grad.lambda  = grad.lambda  + wt * g.lambda;

            mean_muA = mean_muA + wt * muAi;
            mean_muB = mean_muB + wt * muBj;
        end
    end

    result.acceptFraction      = 2 * acos(1 - lambda) / pi;
    result.P_only_D0_avg_valid = result.P_only_D0;
    result.P_only_D1_avg_valid = result.P_only_D1;
    result.P_D0_avg_valid      = result.P_only_D0 + result.P_double;
    result.P_D1_avg_valid      = result.P_only_D1 + result.P_double;

    result.lambda = lambda;
    result.tauA   = tauA;
    result.tauB   = tauB;
    result.etaA   = tauA;
    result.etaB   = tauB;

    result.mean_muA = mean_muA;
    result.mean_muB = mean_muB;
    result.mean_a   = mean_muA * tauA;
    result.mean_b   = mean_muB * tauB;

    result.eventConvention = ...
        'Fluctuation-averaged XX-screened: P_only_D0 = Q10, P_only_D1 = Q01';
end


% =========================================================================
function result = init_result_XX_fluct()
    result = struct();
    result.P_only_D0 = 0;
    result.P_only_D1 = 0;
    result.P_none    = 0;
    result.P_double  = 0;
    result.P_rest    = 0;

    result.P_only_D0_valid = 0;
    result.P_only_D1_valid = 0;
    result.P_none_valid    = 0;
    result.P_double_valid  = 0;
    result.P_rest_valid    = 0;
end


% =========================================================================
function g = init_grad_rowvec_fluct_XX()
    g.muA         = [0, 0, 0];
    g.muB         = [0, 0, 0];
    g.fluct_sigma = [0, 0, 0];
    g.lambda      = [0, 0, 0];

    g.pd0     = [0, 0, 0];
    g.pc0     = [0, 0, 0];
    g.pa0     = [0, 0, 0];
    g.pd1     = [0, 0, 0];
    g.pc1     = [0, 0, 0];
    g.pa1     = [0, 0, 0];

    g.alpha_A = [0, 0, 0];
    g.alpha_B = [0, 0, 0];
    g.dAC     = [0, 0, 0];
    g.dBC     = [0, 0, 0];
end