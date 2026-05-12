function thetaCell = theta2cell(thetaR, varR, colSizes)
% theta2cell: Convert parameter matrix to cell array based on variable names and column sizes.
%
% Inputs:
%   thetaR    - Matrix of parameter values (n rows, m columns). Each column corresponds
%               to a specific variable instance (e.g., different intensity values).
%   varR      - Cell array of variable names (1 x nVars). The order must match the
%               original grouping of columns in thetaR.
%   colSizes  - (Optional) Vector of positive integers of length nVars, specifying
%               how many consecutive columns in thetaR belong to each variable.
%               If not provided, it is assumed that each variable occupies exactly one column,
%               i.e., colSizes = ones(1, length(varR)).
%
% Output:
%   thetaCell - Cell array of length nVars. Each cell contains an n x colSizes(i) matrix
%               of values for the corresponding variable.
%
% Example:
%   thetaR = [0.1 0.2 0.3 0.4 0.5 0.6; 
%             0.7 0.8 0.9 1.0 1.1 1.2];  % 2 rows, 6 columns
%   varR = {'signal_A', 'decoy_A', 'alpha_A', 'dAC', 'signal_B', 'decoy_B', 'alpha_B', 'dBC'};
%   colSizes = [2, 1, 1, 1, 2, 1, 1, 1];  % signal_A and signal_B occupy 2 columns each
%   thetaCell = theta2cell(thetaR, varR, colSizes);
%   thetaCell{1}  % returns 2x2 matrix for signal_A
%
% Copyright (c) 2024 YourName
% Licensed under the MIT License (see LICENSE file for full details).

    % Get number of rows
    n = size(thetaR, 1);

    % Number of variables
    nVars = length(varR);

    % If colSizes not provided, assume each variable occupies one column
    if nargin < 3
        colSizes = ones(1, nVars);
    end

    % Validate total columns
    if sum(colSizes) ~= size(thetaR, 2)
        error('Sum of colSizes (%d) does not match number of columns in thetaR (%d).', ...
              sum(colSizes), size(thetaR, 2));
    end

    % Split thetaR by column blocks
    thetaCell = mat2cell(thetaR, n, colSizes);

end