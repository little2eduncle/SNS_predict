function outcome = SNS_protocol_core(N, thetaA, thetaB, thetaC, opt)

% --- 参数 ---
signal_A = thetaA{1};
decoy_A  = thetaA{2};
alpha_A  = thetaA{3};
dAC      = thetaA{4};

signal_B = thetaB{1};
decoy_B  = thetaB{2};
alpha_B  = thetaB{3};
dBC      = thetaB{4};

pa0 = thetaC{1}; pa1 = thetaC{2};
pc0 = thetaC{3}; pc1 = thetaC{4};
pd0 = thetaC{5}; pd1 = thetaC{6};

% 信道
etaA = 10^(-alpha_A * dAC / 10);
etaB = 10^(-alpha_B * dBC / 10);

% 输出
outcome = zeros(1, N);

for i = 1:N
    
    % === 1. 窗口选择 ===
    A_isZ = rand < opt.z_window;
    B_isZ = rand < opt.z_window;
    
    % === 2. 状态选择 ===
    [muA, catA] = sample_state(A_isZ, signal_A, decoy_A, opt);
    [muB, catB] = sample_state(B_isZ, signal_B, decoy_B, opt);
    
    % === 3. 相位 ===
    delta = 2*pi*rand;
    
    % === 4. 波动 ===
    muA = fluct(muA, opt.fluct_sigma);
    muB = fluct(muB, opt.fluct_sigma);
    
    % === 5. 信道 ===
    muA = muA * etaA * pc0;
    muB = muB * etaB * pc0;
    
    % === 6. 干涉 ===
    I0 = 0.5*(muA + muB + 2*sqrt(muA*muB)*cos(delta));
    I1 = 0.5*(muA + muB - 2*sqrt(muA*muB)*cos(delta));
    
    % === 7. 点击概率 ===
    p0 = 1 - (1-pd0)*(1-pa0)*exp(-I0);
    p1 = 1 - (1-pd1)*(1-pa1)*exp(-I1);
    
    % === 8. 采样 ===
    D0 = rand < p0;
    D1 = rand < p1;
    
    outcome(i) = 1 + D1 + 2*D0;
end

end