function  Pout = fiberlaser_twoend
    global R1 R2 Ppl Ppr sigma_ap sigma_ep sigma_as sigma_es qe N Pssat Ppsat L mu h nu_s 
    
    %----------------参数设置----------------
    lambda_s = 1550 * 1e-9; %
    lambda_p = 976e-9; % 泵浦波长（m）
    tau = 8e-3; % 激发离子的自发衰变速率(亚稳态寿命)
    A_c = 3.1416e-10; % 纤芯截面积
    sigma_ap = 1.7e-25; % 泵浦吸收面积（m^2）
    sigma_ep = 0; % 泵浦发射面积（m^2）
    [sigma_as,sigma_es] = GetErSpectrum(lambda_s);
    N = 6e25; % 单位体积掺杂离子数（m^-3）
    qe = 0.85; % 泵浦量子效率
    R1 =.99;
    R2 =.035;
    L = 40;
    mu = 2; % 偏振模数目
    
    %----------------物理常数及中间过程参数----------------
    c = 3e8; % 真空中光速
    h = 6.626e-34; % 普朗克常量
    nu_s = c/lambda_s; % 
    nu_p = c/lambda_p; % 泵浦光频率   
    Pth = h * nu_s * A_c * N / (tau * pi); % 阈值功率
    Pssat = h * nu_s * A_c * N /( qe * (sigma_es+sigma_as) * tau); % 饱和功率
    Ppsat = h * nu_p * A_c * N /( qe * (sigma_ep+sigma_ap) * tau);
    
    %抽运光功率设置
    Ppl = 100;
    Ppr = 0;
    
    %----------------边值问题数值求解----------------
    OPTION = bvpset('Stats','ON');
    solinit = bvpinit(linspace(0,L,100),[Ppl Ppr 30 Ppr]);
    sol = bvp4c(@f,@fsbc,solinit,OPTION);
    
    %----------------数值计算结果分析和显示----------------
    x = [sol.x];
    y = [sol.y];
    Pout = y(3,end)*(1-R2);
   
    figure;
    subplot(2,1,1)
    plot(x,y(1,:),'b.-',x,y(2,:),'g*-',x,y(3,:),'r',x,y(4,:),'k--');
    grid on;
    title('Pump and laser powers');
    legend('Pp+(z)','Pp-(z)','Ps+(z)','Ps-(z)');
    xlabel('Position z (m)');
    ylabel('Power (W)');
    subplot(2,1,2)
    plot(x,nz)
    grid on;
    title('Relative population density')
    xlabel('Position z (m)');
    ylabel('N_2/N');
   
    
%----------------由目标波长计算光纤激光发射/吸收面积----------------    
function [emint,abint] = GetErSpectrum(wlint)
    wlint = wlint*1e9;
    wl = [1460:10:1520, 1525:5:1550, 1560:10:1640];
    ab = 1E-24*[0.025, 0.05, 0.15, 0.27, 0.24, 0.20, 0.21... % 1460 to 1520
        0.25, 0.4, 0.6, 0.39, 0.29, 0.24, ... % 1525 to 1550
        0.14, 0.06, 0.04, 0.03, 0.02, 0.015, 0.01, 0.002, 0.001]; % 1560 to 1640
    em = 1E-24*[0.001, 0.02, 0.04, 0.1, 0.12, 0.14, 0.17... % 1460 to 1520
        0.2, 0.4, 0.6, 0.42, 0.35, 0.3, ... % 1525 to 1550
        0.24, 0.13, 0.09, 0.1, 0.09, 0.08, 0.07, 0.04, 0.02]; % 1560 to 1640
    if nargin<1,
        wlint = 1460:0.5:1640; 
    end
    abint = interp1(wl,ab,wlint,'PCHIP');%分段三次Hermite插值
    emint = interp1(wl,em,wlint,'PCHIP');
    
    
%----------------由目标波长计算光纤激光发射/吸收面积----------------
function dy = f(x,y)
    global sigma_ap sigma_ep sigma_as sigma_es qe N Pssat Ppsat mu h nu_s 
    delta_nu_s = ;
    dy = zeros(3,1); 
    Nz=N*(sigma_ap/(sigma_ap+sigma_ep)*(y(1)+y(2))/Ppsat+sigma_as/...
     (sigma_as+sigma_es)*(y(3)+y(4))/Pssat)/((y(1)+y(2))/Ppsat+1+(y(4)+y(3))/Pssat); % 平均上能级粒子数浓度与铒离子平均总掺杂浓度之比
 
    dy(1) = (qe*(sigma_ap+sigma_ep)*Nz - qe*sigma_ap) * y(1);
    dy(2) = (qe*(sigma_ap+sigma_ep)*Nz - qe*sigma_ap) * y(2) + mu * qe*sigma_ep*Nz*h*nu_s*delta_nu_s;
    dy(3) = qe*sigma_ap*Nz * y(3);

    
%----------------边界条件----------------
function res = fsbc(y0,yL)
    global R1 R2 Ppl Ppr
    res = [y0(1)-Ppl
       yL(2)-Ppr
       y0(3)-R1*y0(4)
       yL(4)-R2*yL(3)];