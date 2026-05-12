function check_gradients_Zwindow_click_prob_asym()
% Finite-difference gradient check for:
%   [result, grad] = Zwindow_click_prob_asym(...)
%
% It compares analytic gradients against finite differences
% for:
%   [P_only_D0, P_only_D1, P_rest]
%
% Works correctly at boundary points such as muB = 0 by
% automatically switching to one-sided differences when needed.

    clc;

    % -------------------------------------------------------------
    % 1) Base point to test
    % -------------------------------------------------------------
    p.muA     = 0.35;
    p.muB     = 0.00;

    p.pd0     = 1e-6;
    p.pc0     = 0.20;
    p.pa0     = 0.01;

    p.pd1     = 1e-6;
    p.pc1     = 0.20;
    p.pa1     = 0.01;

    p.alpha_A = 0.2;
    p.alpha_B = 0.2;
    p.dAC     = 50;
    p.dBC     = 50;

    % -------------------------------------------------------------
    % 2) Parameters to check
    % -------------------------------------------------------------
    param_names = { ...
        'muA', 'muB', ...
        'pd0', 'pc0', 'pa0', ...
        'pd1', 'pc1', 'pa1', ...
        'alpha_A', 'alpha_B', ...
        'dAC', 'dBC'};

    % -------------------------------------------------------------
    % 3) Evaluate analytic result/gradient at base point
    % -------------------------------------------------------------
    [result0, grad0] = call_model(p);

    fprintf('Base probabilities:\n');
    fprintf('  P_only_D0 = %.16g\n', result0.P_only_D0);
    fprintf('  P_only_D1 = %.16g\n', result0.P_only_D1);
    fprintf('  P_rest    = %.16g\n', result0.P_rest);
    fprintf('  Sum check = %.16g\n\n', ...
        result0.P_only_D0 + result0.P_only_D1 + result0.P_rest);

    % -------------------------------------------------------------
    % 4) Finite-difference check
    % -------------------------------------------------------------
    rows = {};
    tol_abs = 1e-6;
    tol_rel = 1e-4;

    for i = 1:numel(param_names)
        name = param_names{i};
        x0   = p.(name);

        h = choose_step(x0);

        [p_plus, p_minus, ok, scheme] = make_fd_points(p, name, h);
        if ~ok
            warning('Skip parameter %s because a valid FD step could not be constructed.', name);
            continue;
        end

        [r_plus,  ~] = call_model(p_plus);
        [r_minus, ~] = call_model(p_minus);

        dx = p_plus.(name) - p_minus.(name);

        fd = zeros(1,3);
        fd(1) = (r_plus.P_only_D0 - r_minus.P_only_D0) / dx;
        fd(2) = (r_plus.P_only_D1 - r_minus.P_only_D1) / dx;
        fd(3) = (r_plus.P_rest    - r_minus.P_rest   ) / dx;

        an = grad0.(name);

        abs_err = abs(an - fd);
        rel_err = abs_err ./ max(1, abs(fd));

        rows(end+1, :) = {name, scheme, ...
            an(1), fd(1), abs_err(1), rel_err(1), ...
            an(2), fd(2), abs_err(2), rel_err(2), ...
            an(3), fd(3), abs_err(3), rel_err(3)}; %#ok<AGROW>
    end

    % -------------------------------------------------------------
    % 5) Print formatted table
    % -------------------------------------------------------------
    fprintf('Gradient check results\n');
    fprintf(['%-10s %-9s | %-12s %-12s %-12s %-12s | ', ...
             '%-12s %-12s %-12s %-12s | ', ...
             '%-12s %-12s %-12s %-12s\n'], ...
             'param', 'scheme', ...
             'an_D0', 'fd_D0', 'absErr', 'relErr', ...
             'an_D1', 'fd_D1', 'absErr', 'relErr', ...
             'an_rest', 'fd_rest', 'absErr', 'relErr');
    fprintf('%s\n', repmat('-',1,172));

    for i = 1:size(rows,1)
        fprintf(['%-10s %-9s | %-12.4e %-12.4e %-12.4e %-12.4e | ', ...
                 '%-12.4e %-12.4e %-12.4e %-12.4e | ', ...
                 '%-12.4e %-12.4e %-12.4e %-12.4e\n'], ...
                 rows{i,1}, rows{i,2}, ...
                 rows{i,3}, rows{i,4}, rows{i,5}, rows{i,6}, ...
                 rows{i,7}, rows{i,8}, rows{i,9}, rows{i,10}, ...
                 rows{i,11}, rows{i,12}, rows{i,13}, rows{i,14});
    end

    fprintf('\nPass/fail summary\n');
    for i = 1:size(rows,1)
        name = rows{i,1};

        abs_ok = all([rows{i,5}, rows{i,9}, rows{i,13}] < tol_abs);
        rel_ok = all([rows{i,6}, rows{i,10}, rows{i,14}] < tol_rel);

        if abs_ok || rel_ok
            status = 'OK';
        else
            status = 'CHECK';
        end

        fprintf('  %-10s : %s\n', name, status);
    end
end


% =====================================================================
function [result, grad] = call_model(p)
    [result, grad] = Zwindow_click_prob_asym( ...
        p.muA, p.muB, ...
        p.pd0, p.pc0, p.pa0, ...
        p.pd1, p.pc1, p.pa1, ...
        p.alpha_A, p.alpha_B, p.dAC, p.dBC);
end


% =====================================================================
function h = choose_step(x)
% Relative step with a minimum floor
    h = 1e-7 * max(1, abs(x));
end


% =====================================================================
function [p_plus, p_minus, ok, scheme] = make_fd_points(p, name, h)
% Build finite-difference points while respecting parameter domains.
% Returns either:
%   scheme = 'central'
% or
%   scheme = 'forward'

    p_plus  = p;
    p_minus = p;
    ok = true;
    scheme = 'central';

    x0 = p.(name);

    switch name
        case {'muA','muB','pc0','pc1','alpha_A','alpha_B','dAC','dBC'}
            % domain: [0, +inf)
            if x0 - h < 0
                % fall back to forward difference at the boundary
                p_plus.(name)  = x0 + h;
                p_minus.(name) = x0;
                scheme = 'forward';
            else
                p_plus.(name)  = x0 + h;
                p_minus.(name) = x0 - h;
            end

        case {'pd0','pd1','pa0','pa1'}
            % domain: [0,1]
            if x0 - h < 0
                p_plus.(name)  = min(1, x0 + h);
                p_minus.(name) = x0;
                scheme = 'forward';
                if p_plus.(name) == p_minus.(name)
                    ok = false;
                end
            elseif x0 + h > 1
                p_plus.(name)  = x0;
                p_minus.(name) = max(0, x0 - h);
                scheme = 'forward';
                if p_plus.(name) == p_minus.(name)
                    ok = false;
                end
            else
                p_plus.(name)  = x0 + h;
                p_minus.(name) = x0 - h;
            end

        otherwise
            p_plus.(name)  = x0 + h;
            p_minus.(name) = x0 - h;
    end
end