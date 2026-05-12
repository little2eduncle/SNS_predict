function [result, grad] = Zwindow_click_prob_asym( ...
    muA, muB, ...
    pd0, pc0, pa0, ...
    pd1, pc1, pa1, ...
    alpha_A, alpha_B, dAC, dBC)
% Zwindow_click_prob_asym
%
% Paper-style analytic click probabilities for a 2-detector SNS/TF-QKD
% model with:
%   - ideal 50:50 BS
%   - full phase averaging: phi ~ Uniform(0, 2*pi)
%   - constant afterpulse model through pa0, pa1
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
% Averaged probabilities:
%   Q00 = P_none
%   Q10 = P_only_D0
%   Q01 = P_only_D1
%   Q11 = P_double
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
    % Optical intensities arriving at Charlie before detector efficiency
    % ---------------------------------------------------------------------
    a = muA * tauA;
    b = muB * tauB;

    da = init_grad_scalar_asym();
    db = init_grad_scalar_asym();

    da.muA     = tauA;
    da.alpha_A = muA * d_tauA_d_alphaA;
    da.dAC     = muA * d_tauA_d_dAC;

    db.muB     = tauB;
    db.alpha_B = muB * d_tauB_d_alphaB;
    db.dBC     = muB * d_tauB_d_dBC;

    % A = (a+b)/2
    A = 0.5 * (a + b);
    dA = init_grad_scalar_asym();

    fields = fieldnames(dA);
    for k = 1:numel(fields)
        f = fields{k};
        dA.(f) = 0.5 * (da.(f) + db.(f));
    end

    % s = a*b, B = sqrt(s)
    s = a * b;
    B = sqrt(max(s, 0));

    ds_grad = init_grad_scalar_asym();
    for k = 1:numel(fields)
        f = fields{k};
        ds_grad.(f) = b * da.(f) + a * db.(f);
    end

    % ---------------------------------------------------------------------
    % Detector-specific total efficiencies and paper-style parameters
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

    dgamma0 = init_grad_scalar_asym();
    dgamma1 = init_grad_scalar_asym();

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

    dC0 = init_grad_scalar_asym();
    dC1 = init_grad_scalar_asym();

    dC0.pd0 = -(1 - pa0);
    dC0.pa0 = -(1 - pd0);

    dC1.pd1 = -(1 - pa1);
    dC1.pa1 = -(1 - pd1);

    % ---------------------------------------------------------------------
    % Stable kernels:
    %   K0 = exp(-gamma0) * I0(beta0)
    %   K1 = exp(-gamma1) * I0(beta1)
    %   Kd = exp(-(gamma0+gamma1)) * I0(beta0-beta1)
    % ---------------------------------------------------------------------
    [K0, G0] = stable_exp_neg_g_I0(gamma0, beta0);
    [K1, G1] = stable_exp_neg_g_I0(gamma1, beta1);
    [Kd, Gd] = stable_exp_neg_g_I0(gamma0 + gamma1, beta0 - beta1);

    T0 = C0 * K0;
    T1 = C1 * K1;
    Td = C0 * C1 * Kd;

    % ---------------------------------------------------------------------
    % Averaged probabilities
    % ---------------------------------------------------------------------
    P_none    = Td;
    P_only_D0 = T1 - Td;
    P_only_D1 = T0 - Td;
    P_double  = 1 - P_none - P_only_D0 - P_only_D1;

    % clip tiny numerical negatives/overshoots
    P_none    = max(min(P_none,    1), 0);
    P_only_D0 = max(min(P_only_D0, 1), 0);
    P_only_D1 = max(min(P_only_D1, 1), 0);
    P_double  = max(min(P_double,  1), 0);

    % renormalize
    Sprob = P_none + P_only_D0 + P_only_D1 + P_double;
    if Sprob > 0
        P_none    = P_none    / Sprob;
        P_only_D0 = P_only_D0 / Sprob;
        P_only_D1 = P_only_D1 / Sprob;
        P_double  = P_double  / Sprob;
    else
        P_none    = 1;
        P_only_D0 = 0;
        P_only_D1 = 0;
        P_double  = 0;
    end

    P_rest = P_none + P_double;
    P_D0   = P_only_D0 + P_double;
    P_D1   = P_only_D1 + P_double;

    % ---------------------------------------------------------------------
    % Gradients
    %
    % Note:
    % Returned gradients correspond to the raw formulas before clipping /
    % renormalization. They are reliable away from clipping-active regions.
    % ---------------------------------------------------------------------
    grad = init_grad_rowvec_asym();

    for k = 1:numel(fields)
        f = fields{k};

        % ---- dK0 = d[exp(-gamma0) I0(pc0*sqrt(s))]
        dc0 = 0;
        if strcmp(f, 'pc0')
            dc0 = 1;
        end
        dI0_beta0 = dI0_csqrt_stable(pc0, dc0, s, ds_grad.(f));
        dK0 = -K0 * dgamma0.(f) + exp(-gamma0) * dI0_beta0;

        % ---- dK1 = d[exp(-gamma1) I0(pc1*sqrt(s))]
        dc1 = 0;
        if strcmp(f, 'pc1')
            dc1 = 1;
        end
        dI0_beta1 = dI0_csqrt_stable(pc1, dc1, s, ds_grad.(f));
        dK1 = -K1 * dgamma1.(f) + exp(-gamma1) * dI0_beta1;

        % ---- dKd = d[exp(-(gamma0+gamma1)) I0((pc0-pc1)*sqrt(s))]
        cd = pc0 - pc1;
        dcd = 0;
        if strcmp(f, 'pc0')
            dcd = 1;
        elseif strcmp(f, 'pc1')
            dcd = -1;
        end
        dI0_betad = dI0_csqrt_stable(cd, dcd, s, ds_grad.(f));
        dK_d = -Kd * (dgamma0.(f) + dgamma1.(f)) + exp(-(gamma0 + gamma1)) * dI0_betad;

        dT0 = dC0.(f) * K0 + C0 * dK0;
        dT1 = dC1.(f) * K1 + C1 * dK1;
        dTd = (dC0.(f) * C1 + C0 * dC1.(f)) * Kd + C0 * C1 * dK_d;

        dP_only_D0 = dT1 - dTd;
        dP_only_D1 = dT0 - dTd;
        dP_rest    = -dP_only_D0 - dP_only_D1;

        grad.(f) = [dP_only_D0, dP_only_D1, dP_rest];
    end

    % ---------------------------------------------------------------------
    % Output struct
    % ---------------------------------------------------------------------
    result = struct();

    result.P_only_D0 = P_only_D0;
    result.P_only_D1 = P_only_D1;
    result.P_none    = P_none;
    result.P_double  = P_double;
    result.P_rest    = P_rest;

    result.P_only_D0_avg_valid = P_only_D0;
    result.P_only_D1_avg_valid = P_only_D1;
    result.P_D0_avg_valid      = P_D0;
    result.P_D1_avg_valid      = P_D1;

    result.acceptFraction = 1;

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

    result.eventConvention = 'P_only_D0 = Q10 (only D0); P_only_D1 = Q01 (only D1)';

    % legacy-compatible aliases
    result.a  = a;
    result.b  = b;
    result.A  = A;
    result.B  = B;

    result.q0 = C0 * exp(-gamma0);
    result.q1 = C1 * exp(-gamma1);

    result.x0 = beta0;
    result.x1 = beta1;

    result.G0 = G0;
    result.G1 = G1;
    result.Gd = Gd;
end


% =========================================================================
function [K, G] = stable_exp_neg_g_I0(g, x)
% Computes:
%   K = exp(-g) * I0(x)
%   G = I0(x)
% using scaled besseli for stability.

    I0_scaled = besseli(0, x, 1);          % exp(-abs(x)) * I0(x)
    G = exp(abs(x)) * I0_scaled;
    K = exp(-g + abs(x)) * I0_scaled;
end


% =========================================================================
function dI0 = dI0_csqrt_stable(c, dc, s, ds)
% Stable derivative of I0(c*sqrt(s))
%
% x = c*sqrt(s)
%
% d/dθ I0(x)
% = I1(x) * ( sqrt(s)*dc + c/(2*sqrt(s))*ds )
%
% Rewrite to avoid singularity at s=0:
% = I1(x)*B*dc + (c^2/2)*(I1(x)/x)*ds
%
% with I1(x)/x -> 1/2 as x -> 0.

    B = sqrt(max(s, 0));
    x = c * B;

    if abs(x) < 1e-12
        I1x = 0;      % I1(0)=0
        H = 0.5;      % limit I1(x)/x -> 1/2
    else
        I1x = besseli(1, x);
        H = I1x / x;
    end

    dI0 = I1x * B * dc + 0.5 * c^2 * H * ds;
end


% =========================================================================
function g = init_grad_scalar_asym()
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
end


% =========================================================================
function g = init_grad_rowvec_asym()
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
end