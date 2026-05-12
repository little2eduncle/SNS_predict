function [result, grad] = XXwindow_click_prob_asym_screened( ...
    muA, muB, ...
    pd0, pc0, pa0, ...
    pd1, pc1, pa1, ...
    alpha_A, alpha_B, dAC, dBC, ...
    lambda)
% XXwindow_click_prob_asym_screened
%
% Paper-style XX-window click probabilities for a 2-detector SNS/TF-QKD
% model with:
%   - ideal 50:50 BS
%   - phase screening: |cos(phi)| >= 1 - lambda
%   - constant afterpulse model through pa0, pa1
%
% -------------------------------------------------------------------------
% IMPORTANT
%   Returned probabilities:
%       P_only_D0, P_only_D1, P_none, P_double, P_rest
%   are CONDITIONAL on the phase being accepted.
%
%   Additional fields *_valid are unconditional accepted probabilities over
%   the full phase range, i.e. acceptFraction * conditional probability.
%
% -------------------------------------------------------------------------
% Paper-style notation
%
%   tauA = 10^(-alpha_A * dAC / 10)
%   tauB = 10^(-alpha_B * dBC / 10)
%
%   a = muA * tauA
%   b = muB * tauB
%   A = (a+b)/2
%   B = sqrt(a*b)
%
% For detector D0:
%   gamma0 = pc0 * A
%   beta0  = pc0 * B
%
% For detector D1:
%   gamma1 = pc1 * A
%   beta1  = pc1 * B
%
% Instantaneous intensities:
%   I_D0(phi) = gamma0 - beta0 * cos(phi)
%   I_D1(phi) = gamma1 + beta1 * cos(phi)
%
% Click probabilities:
%   P_D0(phi) = 1 - (1-pd0)(1-pa0)exp(-I_D0(phi))
%   P_D1(phi) = 1 - (1-pd1)(1-pa1)exp(-I_D1(phi))
%
% -------------------------------------------------------------------------
% Outputs
%   result.P_*        : conditional on accepted phase
%   result.P_*_valid  : unconditional accepted probabilities
%   result.acceptFraction
%
%   grad.(field) = [dP_only_D0, dP_only_D1, dP_rest]
%   for all parameters including lambda
%
% -------------------------------------------------------------------------

    % ---------------------------------------------------------------------
    % Input checks
    % ---------------------------------------------------------------------
    assert(isscalar(muA) && isscalar(muB), 'muA and muB must be scalars.');
    assert(muA >= 0 && muB >= 0, 'muA and muB must be nonnegative.');

    assert(isscalar(pd0) && isscalar(pd1), 'pd0 and pd1 must be scalars.');
    assert(pd0 >= 0 && pd0 <= 1 && pd1 >= 0 && pd1 <= 1, ...
        'pd0 and pd1 must be in [0,1].');

    assert(isscalar(pa0) && isscalar(pa1), 'pa0 and pa1 must be scalars.');
    assert(pa0 >= 0 && pa0 <= 1 && pa1 >= 0 && pa1 <= 1, ...
        'pa0 and pa1 must be in [0,1].');

    assert(isscalar(pc0) && isscalar(pc1), 'pc0 and pc1 must be scalars.');
    assert(pc0 >= 0 && pc1 >= 0, 'pc0 and pc1 must be nonnegative.');

    assert(isscalar(alpha_A) && isscalar(alpha_B), ...
        'alpha_A and alpha_B must be scalars.');
    assert(alpha_A >= 0 && alpha_B >= 0, ...
        'alpha_A and alpha_B must be nonnegative.');

    assert(isscalar(dAC) && isscalar(dBC), 'dAC and dBC must be scalars.');
    assert(dAC >= 0 && dBC >= 0, 'dAC and dBC must be nonnegative.');

    assert(isscalar(lambda), 'lambda must be a scalar.');
    assert(lambda > 0 && lambda <= 1, ...
        'lambda must satisfy 0 < lambda <= 1.');

    quadAbsTol = 1e-12;
    quadRelTol = 1e-9;
    ln10 = log(10);

    % ---------------------------------------------------------------------
    % Channel transmittance
    % ---------------------------------------------------------------------
    tauA = 10^(-alpha_A * dAC / 10);
    tauB = 10^(-alpha_B * dBC / 10);

    d_tauA_d_alphaA = tauA * (-(ln10 / 10) * dAC);
    d_tauA_d_dAC    = tauA * (-(ln10 / 10) * alpha_A);

    d_tauB_d_alphaB = tauB * (-(ln10 / 10) * dBC);
    d_tauB_d_dBC    = tauB * (-(ln10 / 10) * alpha_B);

    % ---------------------------------------------------------------------
    % Optical intensities before detector efficiency
    % ---------------------------------------------------------------------
    a = muA * tauA;
    b = muB * tauB;
    A = 0.5 * (a + b);

    da = init_grad_scalar_asym_xx();
    db = init_grad_scalar_asym_xx();

    da.muA     = tauA;
    da.alpha_A = muA * d_tauA_d_alphaA;
    da.dAC     = muA * d_tauA_d_dAC;

    db.muB     = tauB;
    db.alpha_B = muB * d_tauB_d_alphaB;
    db.dBC     = muB * d_tauB_d_dBC;

    dA = init_grad_scalar_asym_xx();
    fields = fieldnames(dA);
    for k = 1:numel(fields)
        f = fields{k};
        dA.(f) = 0.5 * (da.(f) + db.(f));
    end

    s = a * b;
    B = sqrt(max(s, 0));

    ds_grad = init_grad_scalar_asym_xx();
    for k = 1:numel(fields)
        f = fields{k};
        ds_grad.(f) = b * da.(f) + a * db.(f);
    end

    % ---------------------------------------------------------------------
    % Detector-specific paper-style parameters
    % ---------------------------------------------------------------------
    etaA0 = pc0 * tauA;
    etaB0 = pc0 * tauB;
    etaA1 = pc1 * tauA;
    etaB1 = pc1 * tauB;

    gamma_A0 = sqrt(max(muA * etaA0, 0));
    gamma_B0 = sqrt(max(muB * etaB0, 0));
    gamma_A1 = sqrt(max(muA * etaA1, 0));
    gamma_B1 = sqrt(max(muB * etaB1, 0));

    gamma0 = pc0 * A;
    gamma1 = pc1 * A;

    beta0  = pc0 * B;
    beta1  = pc1 * B;

    dgamma0 = init_grad_scalar_asym_xx();
    dgamma1 = init_grad_scalar_asym_xx();

    for k = 1:numel(fields)
        f = fields{k};
        dgamma0.(f) = pc0 * dA.(f);
        dgamma1.(f) = pc1 * dA.(f);
    end
    dgamma0.pc0 = A;
    dgamma1.pc1 = A;

    % ---------------------------------------------------------------------
    % Constant prefactors from dark count and afterpulse
    % ---------------------------------------------------------------------
    C0 = (1 - pd0) * (1 - pa0);
    C1 = (1 - pd1) * (1 - pa1);

    dC0 = init_grad_scalar_asym_xx();
    dC1 = init_grad_scalar_asym_xx();

    dC0.pd0 = -(1 - pa0);
    dC0.pa0 = -(1 - pd0);

    dC1.pd1 = -(1 - pa1);
    dC1.pa1 = -(1 - pd1);

    % ---------------------------------------------------------------------
    % Phase acceptance set
    %   |cos(phi)| >= 1 - lambda
    % ---------------------------------------------------------------------
    cthr = 1 - lambda;
    Delta = acos(cthr);
    dDelta_dlambda = 1 / sqrt(lambda * (2 - lambda));

    acceptFraction = 2 * Delta / pi;

    % ---------------------------------------------------------------------
    % Screened kernels:
    %   M(x,Delta) = conditional average of exp(x cos(phi)) over accepted set
    %   N(x,Delta) = dM/dx
    % ---------------------------------------------------------------------
    [M0, N0, dM0_dDelta, H0] = screened_kernel_MN(beta0, Delta, quadAbsTol, quadRelTol);
    [M1, N1, dM1_dDelta, H1] = screened_kernel_MN(beta1, Delta, quadAbsTol, quadRelTol);
    [Md, Nd, dMd_dDelta, Hd] = screened_kernel_MN(beta0 - beta1, Delta, quadAbsTol, quadRelTol);

    exp_neg_g0  = exp(-gamma0);
    exp_neg_g1  = exp(-gamma1);
    exp_neg_gd  = exp(-(gamma0 + gamma1));

    T0 = C0 * exp_neg_g0 * M0;
    T1 = C1 * exp_neg_g1 * M1;
    Td = C0 * C1 * exp_neg_gd * Md;

    % ---------------------------------------------------------------------
    % Conditional accepted probabilities
    % ---------------------------------------------------------------------
    P_none    = Td;
    P_only_D0 = T1 - Td;
    P_only_D1 = T0 - Td;
    P_double  = 1 - P_none - P_only_D0 - P_only_D1;
    P_rest    = P_none + P_double;

    % ---------------------------------------------------------------------
    % Unconditional accepted probabilities over full phase
    % ---------------------------------------------------------------------
    P_none_valid    = acceptFraction * P_none;
    P_only_D0_valid = acceptFraction * P_only_D0;
    P_only_D1_valid = acceptFraction * P_only_D1;
    P_double_valid  = acceptFraction * P_double;
    P_rest_valid    = acceptFraction * P_rest;

    P_D0 = P_only_D0 + P_double;
    P_D1 = P_only_D1 + P_double;

    % ---------------------------------------------------------------------
    % Gradients (semi-analytic)
    % ---------------------------------------------------------------------
    grad = init_grad_rowvec_asym_xx();

    for k = 1:numel(fields)
        f = fields{k};

        if strcmp(f, 'lambda')
            % lambda enters only through Delta
            dM0 = dM0_dDelta * dDelta_dlambda;
            dM1 = dM1_dDelta * dDelta_dlambda;
            dMd = dMd_dDelta * dDelta_dlambda;

            dT0 = C0 * exp_neg_g0 * dM0;
            dT1 = C1 * exp_neg_g1 * dM1;
            dTd = C0 * C1 * exp_neg_gd * dMd;

        else
            % dM(c*sqrt(s),Delta) with fixed Delta
            dc0 = double(strcmp(f, 'pc0'));
            dc1 = double(strcmp(f, 'pc1'));
            dcd = dc0 - dc1;

            dM0 = dM_csqrt_screened(pc0, dc0, s, ds_grad.(f), Delta, N0, H0);
            dM1 = dM_csqrt_screened(pc1, dc1, s, ds_grad.(f), Delta, N1, H1);
            dMd = dM_csqrt_screened(pc0 - pc1, dcd, s, ds_grad.(f), Delta, Nd, Hd);

            dT0 = dC0.(f) * exp_neg_g0 * M0 + C0 * exp_neg_g0 * (dM0 - M0 * dgamma0.(f));
            dT1 = dC1.(f) * exp_neg_g1 * M1 + C1 * exp_neg_g1 * (dM1 - M1 * dgamma1.(f));
            dTd = (dC0.(f) * C1 + C0 * dC1.(f)) * exp_neg_gd * Md ...
                + C0 * C1 * exp_neg_gd * (dMd - Md * (dgamma0.(f) + dgamma1.(f)));
        end

        dP_only_D0 = dT1 - dTd;
        dP_only_D1 = dT0 - dTd;
        dP_rest    = -dP_only_D0 - dP_only_D1;

        grad.(f) = [dP_only_D0, dP_only_D1, dP_rest];
    end

    % ---------------------------------------------------------------------
    % Output
    % ---------------------------------------------------------------------
    result = struct();

    % conditional on accepted phase
    result.P_only_D0 = P_only_D0;
    result.P_only_D1 = P_only_D1;
    result.P_none    = P_none;
    result.P_double  = P_double;
    result.P_rest    = P_rest;

    result.P_only_D0_avg_valid = P_only_D0;
    result.P_only_D1_avg_valid = P_only_D1;
    result.P_D0_avg_valid      = P_D0;
    result.P_D1_avg_valid      = P_D1;

    % unconditional accepted probabilities
    result.P_only_D0_valid = P_only_D0_valid;
    result.P_only_D1_valid = P_only_D1_valid;
    result.P_none_valid    = P_none_valid;
    result.P_double_valid  = P_double_valid;
    result.P_rest_valid    = P_rest_valid;

    result.acceptFraction = acceptFraction;
    result.lambda = lambda;
    result.Delta  = Delta;

    result.tauA = tauA;
    result.tauB = tauB;

    % legacy aliases
    result.etaA = tauA;
    result.etaB = tauB;

    result.etaA0 = etaA0;
    result.etaB0 = etaB0;
    result.etaA1 = etaA1;
    result.etaB1 = etaB1;

    result.gamma_A0 = gamma_A0;
    result.gamma_B0 = gamma_B0;
    result.gamma0   = gamma0;
    result.beta0    = beta0;

    result.gamma_A1 = gamma_A1;
    result.gamma_B1 = gamma_B1;
    result.gamma1   = gamma1;
    result.beta1    = beta1;

    result.M0 = M0;
    result.M1 = M1;
    result.Md = Md;
    result.N0 = N0;
    result.N1 = N1;
    result.Nd = Nd;

    result.eventConvention = ...
        'Conditional on accepted phase: P_only_D0 = Q10, P_only_D1 = Q01';

    % legacy-compatible aliases
    result.a  = a;
    result.b  = b;
    result.A  = A;
    result.B  = B;

    result.x0 = beta0;
    result.x1 = beta1;
end


% =========================================================================
function [M, N, dMdDelta, H0] = screened_kernel_MN(x, Delta, quadAbsTol, quadRelTol)
% M(x,Delta) = conditional average of exp(x cos(phi)) over accepted set
%            = (1/Delta) * integral_0^Delta cosh(x cos(phi)) dphi
%
% N(x,Delta) = dM/dx
%            = (1/Delta) * integral_0^Delta cos(phi) sinh(x cos(phi)) dphi
%
% dM/dDelta  = [cosh(x cos(Delta)) - M] / Delta
%
% H0         = limit of N(x,Delta)/x as x->0
%            = conditional average of cos^2(phi) over accepted set

    funM = @(phi) cosh(x .* cos(phi));
    funN = @(phi) cos(phi) .* sinh(x .* cos(phi));

    IM = integral(funM, 0, Delta, 'ArrayValued', true, ...
        'AbsTol', quadAbsTol, 'RelTol', quadRelTol);
    IN = integral(funN, 0, Delta, 'ArrayValued', true, ...
        'AbsTol', quadAbsTol, 'RelTol', quadRelTol);

    M = IM / Delta;
    N = IN / Delta;

    dMdDelta = (cosh(x * cos(Delta)) - M) / Delta;

    H0 = 0.5 + sin(2 * Delta) / (4 * Delta);
end


% =========================================================================
function dM = dM_csqrt_screened(c, dc, s, ds, Delta, N, H0)
% Stable derivative of M(c*sqrt(s), Delta) wrt a generic parameter.
%
% Let x = c*sqrt(s). Then:
%   dM = N(x,Delta) * sqrt(s) * dc
%      + 0.5 * c^2 * H(x,Delta) * ds
%
% where H(x,Delta) = N(x,Delta) / x, and H(0,Delta) = H0.

    B = sqrt(max(s, 0));
    x = c * B;

    if abs(x) < 1e-12
        H = H0;
    else
        H = N / x;
    end

    dM = N * B * dc + 0.5 * c^2 * H * ds;
end


% =========================================================================
function g = init_grad_scalar_asym_xx()
    g.muA     = 0;
    g.muB     = 0;
    g.pd0     = 0;
    g.pc0     = 0;
    g.pa0     = 0;
    g.pd1     = 0;
    g.pc1     = 0;
    g.pa1     = 0;
    g.alpha_A = 0;
    g.alpha_B = 0;
    g.dAC     = 0;
    g.dBC     = 0;
    g.lambda  = 0;
end


% =========================================================================
function g = init_grad_rowvec_asym_xx()
    g.muA     = [0, 0, 0];
    g.muB     = [0, 0, 0];
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
    g.lambda  = [0, 0, 0];
end