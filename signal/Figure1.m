clear

% Alice和Bob不加载代表基选择及比特选择的相位, 而是以诱骗或信号的选择代表基选择, 以发送或不发送的选择代表比特选择
% 要求 Alice 和 Bob 根据通信双方的诱骗或信号选择对测量结果进行分类

mu=0.35;%真空态平均光子数μ
mv=0.12;%诱饵态平均光子数ν
mvv=0.01;%诱饵态平均光子数ν

%%%%%%%%%%%%%%%%%%%%%%%%%理想光源下的SNS协议

%理想情况下的参数
e_det=0.01; %e_det
p_dark=10^(-11);%单个探测器的暗计数率
yita_D=0.8;%探测器效率
alpha=0.2;%光纤衰减系数α
Y0=p_dark;

%实际情况下的参数
% e_det=0.03; %e_det 光路误码率
% p_dark=10^(-8);%单个探测器的暗计数率
% yita_D=0.3;%探测器效率
% alpha=0.2;%光纤衰减系数α
% Y0=p_dark;% 零光子探测率

N_m=1*10^9;

yita_BS=0.0909;
yita_pnr=0.04127;

len=2000;
Rmax=zeros(1,len);
for k=1:len   %扫描传输距离l

l=500/len*k;
l_a=l;
l_b=l;
         

e0=0.5; %真空态误码率
Q0=Y0;

sigma=0*0.01;

%%%%%%%%%%%%%%%%%%%%%--------------------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 联合态光子数分布
% p0_s_U=p0_s_A_U*p0_s_B_U; p0_s_L=p0_s_A_L*p0_s_B_L;
% p0_d_U=p0_d_A_U*p0_d_B_U; p0_d_L=p0_d_A_L*p0_d_B_L;
% p0_dd_U=p0_dd_A_U*p0_dd_B_U; p0_dd_L=p0_dd_A_L*p0_dd_B_L;
% 
% p1_s_U=p1_s_A_U*p0_s_B_U+p1_s_B_U*p0_s_A_U; p1_s_L=p1_s_A_L*p0_s_B_L+p1_s_B_L*p0_s_A_L; 
% p1_d_U=p1_d_A_U*p0_d_B_U+p1_d_B_U*p0_d_A_U; p1_d_L=p1_d_A_L*p0_d_B_L+p1_d_B_L*p0_d_A_L; 
% p1_dd_U=p1_dd_A_U*p0_dd_B_U+p1_dd_B_U*p0_dd_A_U; p1_dd_L=p1_dd_A_L*p0_dd_B_L+p1_dd_B_L*p0_dd_A_L; 
% 
% p2_s_U=p0_s_A_U*p2_s_B_U+p1_s_A_U*p1_s_B_U+p2_s_A_U*p0_s_B_U; p2_s_L=p0_s_A_L*p2_s_B_L+p1_s_A_L*p1_s_B_L+p2_s_A_L*p0_s_B_L;
% p2_d_U=p0_d_A_U*p2_d_B_U+p1_d_A_U*p1_d_B_U+p2_d_A_U*p0_d_B_U; p2_d_L=p0_d_A_L*p2_d_B_L+p1_d_A_L*p1_d_B_L+p2_d_A_L*p0_d_B_L;
% p2_dd_U=p0_dd_A_U*p2_dd_B_U+p1_dd_A_U*p1_dd_B_U+p2_dd_A_U*p0_dd_B_U; p2_dd_L=p0_dd_A_L*p2_dd_B_L+p1_dd_A_L*p1_dd_B_L+p2_dd_A_L*p0_dd_B_L;

p0_s_U=exp(-2*mu);p0_s_L=exp(-2*mu);
p0_d_U=exp(-2*mv);p0_d_L=exp(-2*mv);
p0_dd_U=exp(-2*mvv);p0_dd_L=exp(-2*mvv);

p1_s_U=2*mu*exp(-2*mu); p1_s_L=2*mu*exp(-2*mu);
p1_d_U=2*mv*exp(-2*mv); p1_d_L=2*mv*exp(-2*mv);
p1_dd_U=2*mvv*exp(-2*mvv); p1_dd_L=2*mvv*exp(-2*mvv);

p2_s_U=2*mu^2*exp(-2*mu); p2_s_L=2*mu^2*exp(-2*mu);
p2_d_U=2*mv^2*exp(-2*mv); p2_d_L=2*mv^2*exp(-2*mv);
p2_dd_U=2*mvv^2*exp(-2*mvv); p2_dd_L=2*mvv^2*exp(-2*mvv);

yita=yita_D*10^(-alpha*l/10);
E_M=0;
S_s=1-(1-Y0)^2*exp(-2*yita*mu); E_s=1/2-1/(2*S_s)*(1-Y0)*(exp(-2*mu*yita*(e_det+E_M))-exp(-2*mu*yita*(1-e_det-E_M)));
S_hs=1-(1-Y0)^2*exp(-yita*mu); E_hs=1/2-1/(2*S_hs)*(1-Y0)*(exp(-mu*yita*(e_det+E_M))-exp(-mu*yita*(1-e_det-E_M)));
S_d=1-(1-Y0)^2*exp(-2*yita*mv); E_d=1/2-1/(2*S_d)*(1-Y0)*(exp(-2*mv*yita*(e_det+E_M))-exp(-2*mv*yita*(1-e_det-E_M)));
S_dd=1-(1-Y0)^2*exp(-2*yita*mvv); E_dd=1/2-1/(2*S_dd)*(1-Y0)*(exp(-2*mvv*yita*(e_det+E_M))-exp(-2*mvv*yita*(1-e_det-E_M)));
S_0=1-(1-Y0)^2; E_0=1/2;

S_1_x_L=(p2_d_L*(S_dd-p0_dd_U*S_0)-p2_dd_U*(S_d-p0_d_L*S_0))/(p2_d_L*p1_dd_L-p2_dd_U*p1_d_U);
e_1_x_U=(S_dd*E_dd-p0_dd_L*S_0/2)/(p1_dd_L*S_1_x_L);

%fe=1.1;
fe=1.16;

for i=1:10
    eps=0.01*i;
    S_Z_SS=eps^2*S_s; SE_Z_SS=eps^2*S_s;
    S_Z_NSNS=(1-eps)^2*S_0; SE_Z_NSNS=(1-eps)^2*S_0;
    S_Z_SNS=2*eps*(1-eps)*S_hs; SE_Z_SNS=0;

    S_Z=S_Z_SS+S_Z_NSNS+S_Z_SNS; E_Z=(SE_Z_SS+SE_Z_NSNS+SE_Z_SNS)/S_Z;
    R(k)=max(2*eps*(1-eps)*p1_s_L*S_1_x_L*(1-h(e_1_x_U))-S_Z*fe*h(E_Z), 0);
    if(R(k)>10^(-12))
        l_long=l;
    end
    if(R(k)>Rmax(k))
        Rmax(k)=R(k);
    end
end

end
% F=quad(@(x)e_mu(x,mu,yita,(yita_s*0.091),sigma),0,2*mu);
%E_s=min(F,0.5);
%e1_s_U=min((F/delta1_s_L),0.5);

t=1:len;%画图
ll=1000/len*t;
p1 = plot(ll,R,'r', 'LineWidth', 2);
set(gca,'YScale','log');
hold on