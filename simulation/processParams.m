function [thetaP, thetaR, thetaF, varR, varF] = processParams(thetas, varF, varR, fluct_sigma)
    %processParams: This function processes the input parameters, separating them into 
    %                fixed and random variables. For each random variable, it computes 
    %                the corresponding prior distribution parameters (alphas, betas, 
    %                upper and lower bounds) based on the distribution type (Gamma or Beta).
    %
    % Inputs:
    %     thetas      - Cell array containing the grouped variables: thetaA, thetaB, thetaC
    %     varF        - Cell array of fixed variable names
    %     varR        - Cell array of random variable names
    %     sigma       - Standard deviation for the prior distributions (default: 0)
    %
    % Outputs:
    %     thetaP      - Prior parameters for the random variables (alphas, betas, ub, lb)
    %     thetaR      - True values of the random variables
    %     thetaF      - Values of the fixed variables
    %     varR        - Updated random variable names
    %     varF        - Updated fixed variable names
    %
    % Copyright (c) 2024 Ibrahim Almosallam <ibrahim@almosallam.org>
    % Licensed under the MIT License (see LICENSE file for full details).

    if nargin < 4 || isempty(fluct_sigma)
        fluct_sigma = 0;
    end
    
    mus = struct();  % Structure to store the means of variables

    % Unpack theta parameters into respective groups (thetaA, thetaB, thetaE)
    [thetaA, thetaB, thetaC] = deal(thetas{:});
    [mus.signal_A, mus.decoy_A, mus.alpha_A, mus.dAC] = deal(thetaA{:});
    [mus.signal_B, mus.decoy_B, mus.alpha_B, mus.dBC] = deal(thetaB{:});

    [mus.pa0, mus.pa1, mus.pc0, mus.pc1, mus.pd0, mus.pd1, mus.phi_center] = deal(thetaC{:});

    % Map variable names to distribution types
    distType = containers.Map();
    gammaVars = {'dAC', 'dBC'};
    for i = 1:numel(gammaVars)
        distType(gammaVars{i}) = 'gamma';
    end
    betaVars = {'signal_A', 'decoy_A', 'alpha_A', 'signal_B', 'decoy_B', 'alpha_B', 'pa0', 'pa1', 'pc0', 'pc1', 'pd0', 'pd1'};
    for i = 1:numel(betaVars)
        distType(betaVars{i}) = 'beta';
    end

    % Collect fixed parameters
    thetaF = cell(1, numel(varF));
    for i = 1:numel(varF)
        varName = varF{i};
        thetaF{i} = mus.(varName);
    end

    % Initialize prior parameter arrays
    thetaR = cell(1, numel(varR));  % True values of random variables
    alphas = cell(1, numel(varR));  % Alpha parameters for priors
    betas = cell(1, numel(varR));   % Beta parameters for priors
    ub = cell(1, numel(varR));      % Upper bounds for priors
    lb = cell(1, numel(varR));      % Lower bounds for priors

    % Process random variables
    for i = 1:numel(varR)
        varName = varR{i};
        thetaR{i} = mus.(varName);  % Assign ground truth values

        mu = mus.(varName);  % Mean value for the variable

        % Assign appropriate prior parameters based on distribution type
        switch distType(varName)
            case 'gamma'
                [alpha_param, beta_param, ub_param, lb_param] = gamma_parameters(mu, (mu * fluct_sigma).^2);
            case 'beta'
                [alpha_param, beta_param, ub_param, lb_param] = beta_parameters(mu, (mu * fluct_sigma).^2);
        end

        % Store prior parameters
        alphas{i} = alpha_param;
        betas{i} = beta_param;
        ub{i} = ub_param;
        lb{i} = lb_param;
    end

    % Pack prior parameters into a cell array
    thetaP = {alphas, betas, ub, lb};
end