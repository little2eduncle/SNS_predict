function [thetaR_MAP, info] = MAP_XX(C, varR, thetaF, thetaP, colSizes, varargin)

    try
        opt_options;
    catch
        options = optimoptions('fminunc', ...
            'SpecifyObjectiveGradient', true, ...
            'Display', 'iter');
    end

    if ~exist('options', 'var')
        options = optimoptions('fminunc', ...
            'SpecifyObjectiveGradient', true, ...
            'Display', 'iter');
    end

    p = inputParser;
    addParameter(p, 'lambda', 0.2);
    addParameter(p, 'fluct_sigma', 0.01);
    addParameter(p, 'algorithm', 'quasi-newton');
    addParameter(p, 'maxIters', 1000);

    % multi-start
    addParameter(p, 'nStarts', 10);
    addParameter(p, 'relNoise', 0.05);
    addParameter(p, 'doDiagnostics', true);

    % 新增：thetaR 每个 cell 对应的变量名
    % 例如 {'decoy_A','decoy_B'}
    addParameter(p, 'thetaRNames', {});

    parse(p, varargin{:});

    lambda        = p.Results.lambda;
    fluct_sigma   = p.Results.fluct_sigma;
    algorithm     = p.Results.algorithm;
    maxIters      = p.Results.maxIters;
    nStarts       = p.Results.nStarts;
    relNoise      = p.Results.relNoise;
    doDiagnostics = p.Results.doDiagnostics;
    thetaRNames   = p.Results.thetaRNames;

    [alphas, betas, ub, lb] = deal(thetaP{:});
    dim = numel(alphas);

    % ---------- 初值：先验均值 ----------
    thetaR0_cells = cell(1, dim);
    for i = 1:dim
        if any(isinf(ub{i}))
            thetaR0_cells{i} = alphas{i} ./ betas{i} + lb{i};
        else
            zmean = alphas{i} ./ (alphas{i} + betas{i});
            thetaR0_cells{i} = zmean .* (ub{i} - lb{i}) + lb{i};
        end
    end
    thetaR0 = flatten_cell_row(thetaR0_cells);

    lbv = flatten_cell_row(lb);
    ubv = flatten_cell_row(ub);

    % 如果用户没显式给 thetaRNames，就尝试默认用 varR
    if isempty(thetaRNames)
        if numel(varR) == dim
            thetaRNames = varR;
        else
            error(['thetaRNames is required when numel(varR) ~= numel(thetaR). ', ...
                   'Example: ''thetaRNames'', {''decoy_A'',''decoy_B''}']);
        end
    end

    if numel(thetaRNames) ~= dim
        error('thetaRNames must have the same number of cells as thetaR / thetaP.');
    end

    % ---------- 目标函数 ----------
    logpdf = @(phi) logPDF(phi, thetaF, varR, C, thetaP, -1, colSizes, ...
                           'lambda', lambda,...
                           'fluct_sigma', fluct_sigma);

    options.Algorithm = algorithm;
    options.MaxIterations = maxIters;
    options.MaxFunctionEvaluations = 5 * maxIters;

    % ---------- multi-start ----------
    starts = cell(1, nStarts);
    starts{1} = make_theta_strictly_interior(thetaR0, lbv, ubv);

    for s = 2:nStarts
        theta_try_cells = sampleThetaRParameters(thetaR0_cells, thetaRNames, varR, relNoise, thetaP);
        theta_try = flatten_cell_row(theta_try_cells);

        theta_try = project_theta_row(theta_try, lbv, ubv);
        theta_try = make_theta_strictly_interior(theta_try, lbv, ubv);

        starts{s} = theta_try;
    end

    best.fval = inf;
    best.theta = [];
    best.phi = [];
    best.exitflag = [];
    best.output = [];
    best.grad = [];
    best.hessian = [];

    allRuns = struct('theta0',{},'thetaStar',{},'fval',{}, ...
                     'exitflag',{},'iterations',{},'firstorderopt',{});

    for s = 1:nStarts
        theta_init = starts{s};
        % disp(['Start ', num2str(s), ':']);
        % disp(theta_init);

        phi_init = Phi(theta_init, thetaP);

        [phiStar, fval, exitflag, output, grad, hessian] = fminunc(logpdf, phi_init, options);
        thetaStar = Theta(phiStar, thetaP);

        allRuns(s).theta0 = theta_init;
        allRuns(s).thetaStar = thetaStar;
        allRuns(s).fval = fval;
        allRuns(s).exitflag = exitflag;
        allRuns(s).iterations = output.iterations;
        allRuns(s).firstorderopt = output.firstorderopt;

        if fval < best.fval
            best.fval = fval;
            best.theta = thetaStar;
            best.phi = phiStar;
            best.exitflag = exitflag;
            best.output = output;
            best.grad = grad;
            best.hessian = hessian;
        end
    end

    thetaR_MAP = best.theta;

    info = struct();
    info.bestFval = best.fval;
    info.bestPhi = best.phi;
    info.bestTheta = best.theta;
    info.bestExitflag = best.exitflag;
    info.bestOutput = best.output;
    info.bestGrad = best.grad;
    info.bestGradNorm = norm(best.grad);
    info.bestHessian = best.hessian;
    info.allRuns = allRuns;

    if doDiagnostics && ~isempty(best.hessian)
        Hsym = 0.5 * (best.hessian + best.hessian.');
        info.hessianEigVals = eig(Hsym);
        info.hessianCond = cond(Hsym);
        info.hessianSymErr = norm(best.hessian - best.hessian.', 'fro');
    end
end

function theta = project_theta_row(theta, lbv, ubv)
    theta = theta(:).';
    lbv   = lbv(:).';
    ubv   = ubv(:).';

    theta = max(theta, lbv);
    finiteMask = ~isinf(ubv);
    theta(finiteMask) = min(theta(finiteMask), ubv(finiteMask));
end


function theta = make_theta_strictly_interior(theta, lbv, ubv)
    theta = theta(:).';
    lbv   = lbv(:).';
    ubv   = ubv(:).';

    epsB = 1e-10;

    finiteMask = ~isinf(ubv);
    semiMask   = isinf(ubv);

    theta(finiteMask) = max(theta(finiteMask), lbv(finiteMask) + epsB);
    theta(finiteMask) = min(theta(finiteMask), ubv(finiteMask) - epsB);

    theta(semiMask) = max(theta(semiMask), lbv(semiMask) + epsB);
end


function row = flatten_cell_row(C)
    row = [C{:}];
    row = row(:).';
end

function sample_thetaR_cells = sampleThetaRParameters(thetaR_cells, thetaRNames, varR, noise, thetaP)
% Sample thetaR cells for multi-start in MAP
%
% Inputs:
%   thetaR_cells  : cell array, e.g. {decoy_A, decoy_B}
%   thetaRNames   : cell array of names, e.g. {'decoy_A','decoy_B'}
%   varR          : names of random variables
%   noise         : relative noise level
%   thetaP        : {alphas, betas, ub, lb}
%
% Output:
%   sample_thetaR_cells : sampled thetaR cells

    [~, ~, ub, lb] = deal(thetaP{:});

    dim = numel(thetaR_cells);
    sample_thetaR_cells = cell(1, dim);

    for i = 1:dim
        mu  = thetaR_cells{i};
        ubi = ub{i};
        lbi = lb{i};

        sample_thetaR_cells{i} = sampleValues_MAP( ...
            thetaRNames{i}, mu, noise, ubi, lbi, varR);
    end
end


function samples = sampleValues_MAP(variable, mu, noise, ub, lb, varR)
% Helper function for MAP starts
%
% If variable is not in varR, return mu unchanged.
% Otherwise:
%   - Gamma for semi-bounded [lb, inf)
%   - Beta  for bounded [lb, ub]

    if (~ismember(variable, varR) || noise == 0)
        samples = mu;
        return;
    end

    if isinf(ub)
        % -------- Gamma distribution for semi-bounded variables --------
        % sigma is relative to mu-lb
        mu_shift = mu - lb;
        mu_shift = max(mu_shift, 1e-12);

        sigma = abs(mu_shift) .* noise;
        sigma = max(sigma, 1e-12);

        [alpha, beta] = gamma_parameters(mu_shift, sigma.^2);

        % MATLAB: gamrnd(shape, scale), where scale = 1 / rate
        samples = gamrnd(alpha, 1 ./ beta) + lb;

    else
        % -------- Beta distribution for bounded variables --------
        width = ub - lb;
        if any(width <= 0)
            error('For bounded variables, ub must be strictly greater than lb.');
        end

        z = (mu - lb) ./ width;
        z = min(max(z, 1e-8), 1 - 1e-8);

        sigma = abs(z) .* noise;
        sigma = max(sigma, 1e-8);

        % Make sure variance is valid for Beta
        maxVar = z .* (1 - z);
        sigma2 = min(sigma.^2, 0.99 .* maxVar);

        [alpha, beta] = beta_parameters(z, sigma2);

        samples = betarnd(alpha, beta) .* width + lb;
    end
end