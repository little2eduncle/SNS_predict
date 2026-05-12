function [mu_nodes, w, dmu_dmu0, dmu_dsigma] = gaussian_fluct_rule(mu0, sigma, n)
% 修改版：截断高斯（区间 [(1-sigma)mu0, (1+sigma)mu0]）+ 数值积分
% 保持接口完全一致（无需改主函数）

    if sigma <= 0 || mu0 == 0
        mu_nodes = mu0;
        w = 1;
        dmu_dmu0 = 1;
        dmu_dsigma = 0;
        return;
    end

    % ---- 截断区间
    a = (1 - sigma) * mu0;
    b = (1 + sigma) * mu0;

    % ---- 用均匀节点（最小改动，不引入lgwt）
    x = linspace(a, b, n)';
    dx = (b - a) / (n - 1);

    % ---- 高斯分布（相对波动）
    std = sigma * mu0;
    pdf = exp(-(x - mu0).^2 ./ (2 * std^2));

    % ---- 归一化
    Z = sum(pdf) * dx;
    w = pdf * dx / Z;

    mu_nodes = x;

    % ---- 梯度（保持你原框架可用）
    dmu_dmu0 = ones(size(x));      % 主导项（够用了）
    dmu_dsigma = (x - mu0) / sigma; % 近似梯度

end