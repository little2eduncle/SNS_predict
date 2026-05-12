function [C_total, C_windows, C_intensity, C_XX_intensity, D0, D1] = ...
    simulate_SNS_TF_QKD(N, thetaA, thetaB, thetaC, varargin)
% simulate_SNS_TF_QKD
% SNS-TF-QKD simulation adapted to the main script.
%
% Gaussian fluctuation version (Scheme 2):
%   - Only nonzero intensities fluctuate
%   - Vacuum stays exactly 0
%   - Fluctuation model:
%         mu_out ~ N(mu_in, fluct_sigma^2), for mu_in > 0
%
% Inputs:
%   N       : total number of pulses
%   thetaA  : {signal_A, decoy_A, alpha_A, dAC}
%   thetaB  : {signal_B, decoy_B, alpha_B, dBC}
%   thetaC  : {pa0, pa1, pc0, pc1, pd0, pd1}
%
% Optional parameters:
%   'buffer'      : chunk size, default 1e4
%   'z_window'    : probability of choosing Z window
%   'epsilon'     : probability of sending signal in Z window
%   'p_decoy'     : probability vector for X-window intensity selection
%                   over [decoy_1, decoy_2, ..., vacuum]
%   'lambda'      : phase-slice threshold for valid XX events
%   'fluct_sigma' : Gaussian std for nonzero intensities
%   'vacuum'      : vacuum intensity (default 0)
%
% Outputs:
%   C_total        : 1x4 total counts [00,01,10,11]
%   C_windows      : struct with counts in ZZ/ZX/XZ/XX windows
%   C_intensity    : struct with counts for different intensity combinations
%   C_XX_intensity : struct with valid XX response events under ALL
%                    non-vacuum decoy_A × decoy_B combinations
%                    that pass the phase mask
%                    counts size = [n_decoyA, n_decoyB, 4]
%   D0, D1         : 1xN logical arrays recording detector clicks
%
% Event encoding:
%   outcome = 1 -> 00 (no click)
%   outcome = 2 -> 01 (only D1 clicks)
%   outcome = 3 -> 10 (only D0 clicks)
%   outcome = 4 -> 11 (double click)

% -------------------------------------------------------------------------
% Parse input arguments
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'buffer', 1e4);
addParameter(p, 'z_window', 0.5);
addParameter(p, 'epsilon', 0.99);
addParameter(p, 'p_decoy', []);
addParameter(p, 'lambda', 0.2);
addParameter(p, 'fluct_sigma', 0.0);
addParameter(p, 'vacuum', 0);
parse(p, varargin{:});
opt = p.Results;

buffer      = opt.buffer;
z_window    = opt.z_window;
epsilon     = opt.epsilon;
p_decoy     = opt.p_decoy;
lambda      = opt.lambda;
fluct_sigma = opt.fluct_sigma;
vacuum      = opt.vacuum;

% -------------------------------------------------------------------------
% Unpack thetaA / thetaB / thetaC
% -------------------------------------------------------------------------
signal_A = thetaA{1};   % scalar
decoy_A  = thetaA{2};   % e.g. [u_d1, u_d2]
alpha_A  = thetaA{3};
dAC      = thetaA{4};

signal_B = thetaB{1};   % scalar
decoy_B  = thetaB{2};   % e.g. [u_d1, u_d2]
alpha_B  = thetaB{3};
dBC      = thetaB{4};

pa0 = thetaC{1};
pa1 = thetaC{2};
pc0 = thetaC{3};
pc1 = thetaC{4};
pd0 = thetaC{5};
pd1 = thetaC{6};

% Force row vectors
decoy_A = decoy_A(:).';
decoy_B = decoy_B(:).';

% X-window intensity sets = [decoys, vacuum]
x_set_A = [decoy_A, vacuum];
x_set_B = [decoy_B, vacuum];

% Default decoy selection probabilities
if isempty(p_decoy)
    p_decoy = ones(1, numel(x_set_A)) / numel(x_set_A);
end

% -------------------------------------------------------------------------
% Basic consistency checks
% -------------------------------------------------------------------------
assert(isscalar(signal_A), 'signal_A must be a scalar.');
assert(isscalar(signal_B), 'signal_B must be a scalar.');
assert(abs(sum(p_decoy) - 1) < 1e-12, 'The elements of p_decoy must sum to 1.');
assert(numel(x_set_A) == numel(p_decoy), ...
    'Length mismatch: p_decoy must match [decoy_A, vacuum].');
assert(numel(x_set_B) == numel(p_decoy), ...
    'Length mismatch: p_decoy must match [decoy_B, vacuum].');
assert(numel(x_set_A) == numel(x_set_B), ...
    'Current version assumes Alice and Bob use the same number of X-window intensities.');

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
n_signalA = 2;                 % Z window bookkeeping: [signal, vacuum]
n_signalB = 2;
n_xA      = numel(x_set_A);    % X-window intensities: [decoys..., vacuum]
n_xB      = numel(x_set_B);

n_decoyA  = numel(decoy_A);    % non-vacuum decoys only
n_decoyB  = numel(decoy_B);

% Four outcomes: [00, 01, 10, 11]
C_total = zeros(1, 4);

C_ZZ = zeros(1, 4);
C_ZX = zeros(1, 4);
C_XZ = zeros(1, 4);
C_XX = zeros(1, 4);

% Counts for different intensity combinations
% signal_* indices: 1=signal, 2=vacuum (in Z window)
% decoy_* indices : over x_set_* = [decoy_1, decoy_2, ..., vacuum] (in X window)
C_intensity.signal_signal = zeros(n_signalA, n_signalB, 4);
C_intensity.signal_decoy  = zeros(n_signalA, n_xB, 4);
C_intensity.decoy_signal  = zeros(n_xA, n_signalB, 4);
C_intensity.decoy_decoy   = zeros(n_xA, n_xB, 4);

% Counts under valid XX events, grouped by ALL non-vacuum decoy combinations
% A-side index over decoy_A, B-side index over decoy_B
C_XX_intensity.mu_nominal_A = decoy_A(:);     % n_decoyA x 1
C_XX_intensity.mu_nominal_B = decoy_B(:).';   % 1 x n_decoyB

% outcome order: [00,01,10,11]
%               [none, only D1, only D0, double]
C_XX_intensity.counts       = zeros(n_decoyA, n_decoyB, 4);
C_XX_intensity.D0_only      = zeros(n_decoyA, n_decoyB);
C_XX_intensity.D1_only      = zeros(n_decoyA, n_decoyB);
C_XX_intensity.single_click = zeros(n_decoyA, n_decoyB);
C_XX_intensity.double_click = zeros(n_decoyA, n_decoyB);
C_XX_intensity.total_valid  = zeros(n_decoyA, n_decoyB);

D0 = false(1, N);
D1 = false(1, N);

num_chunks = ceil(N / buffer);

etaA = 10^(-alpha_A * dAC / 10);
etaB = 10^(-alpha_B * dBC / 10);

% -------------------------------------------------------------------------
% Main simulation loop
% -------------------------------------------------------------------------
start_up;

for chunk = 1:num_chunks

    first = (chunk - 1) * buffer + 1;
    last  = min(chunk * buffer, N);
    ind   = first:last;
    chunk_size = last - first + 1;

    % Window choice: true = Z, false = X
    A_isZ = rand(1, chunk_size) < z_window;
    B_isZ = rand(1, chunk_size) < z_window;

    % State preparation
    [muA_nominal, sendA, phaseA, idxA, typeA] = prepare_party_states( ...
        A_isZ, signal_A, vacuum, x_set_A, epsilon, p_decoy);

    [muB_nominal, sendB, phaseB, idxB, typeB] = prepare_party_states( ...
        B_isZ, signal_B, vacuum, x_set_B, epsilon, p_decoy);

    % Gaussian intensity fluctuation (Scheme 2):
    % only nonzero intensities fluctuate; vacuum stays 0
    muA = apply_intensity_fluctuation_local(muA_nominal, fluct_sigma);
    muB = apply_intensity_fluctuation_local(muB_nominal, fluct_sigma);

    % Effective intensities after channel loss
    muA_eff = muA .* etaA;
    muB_eff = muB .* etaB;

    % Phase difference
    delta_theta = phaseA - phaseB;

    % Mean photon numbers at the two output ports of the 50:50 BS
    I0 = 0.5 .* (muA_eff + muB_eff + 2 .* sqrt(muA_eff .* muB_eff) .* cos(delta_theta));
    I1 = 0.5 .* (muA_eff + muB_eff - 2 .* sqrt(muA_eff .* muB_eff) .* cos(delta_theta));

    I0 = max(I0, 0);
    I1 = max(I1, 0);

    % Detector click probabilities
    p0_final = 1 - (1 - pd0) .* (1 - pa0) .* exp(-pc0 .* I0);
    p1_final = 1 - (1 - pd1) .* (1 - pa1) .* exp(-pc1 .* I1);

    p0_final = min(max(p0_final, 0), 1);
    p1_final = min(max(p1_final, 0), 1);

    % Independent Bernoulli sampling
    finalD0_chunk = rand(1, chunk_size) < p0_final;
    finalD1_chunk = rand(1, chunk_size) < p1_final;

    D0(ind) = finalD0_chunk;
    D1(ind) = finalD1_chunk;

    % Phase-slice condition for valid XX events
    phase_mask = (1 - abs(cos(delta_theta))) <= lambda;

    % Count events one by one
    for j = 1:chunk_size

        d0j = finalD0_chunk(j);
        d1j = finalD1_chunk(j);

        % Outcome encoding:
        % 1 -> 00
        % 2 -> 01 (only D1 clicks)
        % 3 -> 10 (only D0 clicks)
        % 4 -> 11
        if ~d0j && ~d1j
            outcome = 1;
        elseif ~d0j && d1j
            outcome = 2;
        elseif d0j && ~d1j
            outcome = 3;
        else
            outcome = 4;
        end

        % Total counts
        C_total(outcome) = C_total(outcome) + 1;

        % Window-wise counts
        if A_isZ(j) && B_isZ(j)

            % In ZZ, keep xor(sendA, sendB)
            if xor(sendA(j), sendB(j))
                C_ZZ(outcome) = C_ZZ(outcome) + 1;
            end

        elseif A_isZ(j) && ~B_isZ(j)

            C_ZX(outcome) = C_ZX(outcome) + 1;

        elseif ~A_isZ(j) && B_isZ(j)

            C_XZ(outcome) = C_XZ(outcome) + 1;

        else
            % In XX, keep ALL non-vacuum decoy-decoy combinations
            % that satisfy the phase-slice condition
            valid_xx_combo = strcmp(typeA{j}, 'decoy') && ...
                             strcmp(typeB{j}, 'decoy') && ...
                             (idxA(j) <= n_decoyA) && ...
                             (idxB(j) <= n_decoyB) && ...
                             phase_mask(j);

            if valid_xx_combo
                C_XX(outcome) = C_XX(outcome) + 1;

                ia = idxA(j);   % 1..n_decoyA
                ib = idxB(j);   % 1..n_decoyB

                C_XX_intensity.counts(ia, ib, outcome) = ...
                    C_XX_intensity.counts(ia, ib, outcome) + 1;

                C_XX_intensity.total_valid(ia, ib) = ...
                    C_XX_intensity.total_valid(ia, ib) + 1;

                if outcome == 2
                    % only D1
                    C_XX_intensity.D1_only(ia, ib) = C_XX_intensity.D1_only(ia, ib) + 1;
                    C_XX_intensity.single_click(ia, ib) = C_XX_intensity.single_click(ia, ib) + 1;

                elseif outcome == 3
                    % only D0
                    C_XX_intensity.D0_only(ia, ib) = C_XX_intensity.D0_only(ia, ib) + 1;
                    C_XX_intensity.single_click(ia, ib) = C_XX_intensity.single_click(ia, ib) + 1;

                elseif outcome == 4
                    % double click
                    C_XX_intensity.double_click(ia, ib) = C_XX_intensity.double_click(ia, ib) + 1;
                end
            end
        end

        % Counts by intensity combination
        key = [typeA{j} '_' typeB{j}];

        switch key
            case 'signal_signal'
                C_intensity.signal_signal(idxA(j), idxB(j), outcome) = ...
                    C_intensity.signal_signal(idxA(j), idxB(j), outcome) + 1;

            case 'signal_decoy'
                C_intensity.signal_decoy(idxA(j), idxB(j), outcome) = ...
                    C_intensity.signal_decoy(idxA(j), idxB(j), outcome) + 1;

            case 'decoy_signal'
                C_intensity.decoy_signal(idxA(j), idxB(j), outcome) = ...
                    C_intensity.decoy_signal(idxA(j), idxB(j), outcome) + 1;

            case 'decoy_decoy'
                C_intensity.decoy_decoy(idxA(j), idxB(j), outcome) = ...
                    C_intensity.decoy_decoy(idxA(j), idxB(j), outcome) + 1;
        end
    end
end

C_windows = struct('ZZ', C_ZZ, ...
                   'ZX', C_ZX, ...
                   'XZ', C_XZ, ...
                   'XX', C_XX);

end

% =========================================================================
% Local helper function 1: state preparation
% =========================================================================
function [mu, send, phase, idx, type] = ...
    prepare_party_states(isZ, signal_mu, vacuum_mu, x_set, epsilon, p_decoy)

n = numel(isZ);

mu    = zeros(1, n);
send  = false(1, n);
phase = 2 * pi * rand(1, n);
idx   = zeros(1, n);
type  = cell(1, n);

cdf_decoy = cumsum(p_decoy);

for j = 1:n
    if isZ(j)
        % Z window: signal with probability epsilon, otherwise vacuum
        if rand < epsilon
            mu(j)   = signal_mu;
            send(j) = true;
            idx(j)  = 1;   % 1 = signal
        else
            mu(j)   = vacuum_mu;
            send(j) = false;
            idx(j)  = 2;   % 2 = vacuum
        end
        type{j} = 'signal';
    else
        % X window: choose from x_set = [decoy_1, decoy_2, ..., vacuum]
        r = rand;
        k = find(r <= cdf_decoy, 1, 'first');
        if isempty(k)
            k = numel(x_set);
        end
        mu(j)   = x_set(k);
        send(j) = false;
        idx(j)  = k;
        type{j} = 'decoy';
    end
end

end

% =========================================================================
% Local helper function 2: Gaussian fluctuation (Scheme 2)
% =========================================================================
function mu_out = apply_intensity_fluctuation_local(mu_in, fluct_sigma)
% Scheme 2:
%   For mu_in > 0:
%       mu_out ~ N(mu_in, fluct_sigma^2)
%   For mu_in = 0:
%       mu_out = 0
%
% After sampling, negative values are truncated to 0.

if fluct_sigma <= 0
    mu_out = mu_in;
    return;
end

mu_out = mu_in;

mask = mu_in > 0;
mu_out(mask) = mu_in(mask) + fluct_sigma .* randn(size(mu_in(mask)));
mu_out(mask) = max(mu_out(mask), 0);

end