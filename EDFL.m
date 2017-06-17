function  Pout = fiberlaser_twoend
    global R1 R2 Ppl Ppr sigma_ap sigma_ep sigma_as sigma_es gamma_s ...
        gamma_p N alpha_p alpha_s Pssat Ppsat Ppsp Ppsv mu kappa elta
    
    %----------------参数设置----------------
    lambda_s = 1100 * 1e-9; %
    lambda_p = 976e-9; % 泵浦光波长
    tau = 10e-3; % Er离子亚稳态寿命
    A_c = 3.1416e-10; % 纤芯截面积
    sigma_ap = 26e-21*1e-4;
    sigma_ep = 26e-21*1e-4;
    sigma_as = 1e-23*1e-4;
    sigma_es = 1.6e-21*1e-4;
    N = 5.5351e+025;
    alpha_p = 2e-5*1e2;
    alpha_s = 4e-6*1e2;
    gamma_s = 0.82;
    gamma_p = 0.0024;
    R1 =.99;
    R2 =.035;
    L = 40;
    
    %----------------物理常数及中间过程参数----------------
    c = 3e8; % 真空中光速
    h = 6.626e-34; % 普朗克常量
    nu_s = c/lambda_s; % 
    nu_p = c/lambda_p; % 泵浦光频率   
    Pth = h * nu_s * A_c * N / tau; % 阈值功率
    Pssat = h * nu_s * A_c/( gamma_s * (sigma_es+sigma_as) * tau); % 饱和功率
    Ppsat = h * nu_p * A_c/( gamma_p * (sigma_ep+sigma_ap) * tau);
    
    %抽运光功率设置
    Ppl = 50;
    Ppr = 50;
    
    %`端面抽运的光纤激光器边值问题数值求解
    OPTION = bvpset('Stats','ON');
    solinit = bvpinit(linspace(0,L,10),[Ppl Ppr 30 Ppr]);
    sol = bvp4c(@f,@fsbc,solinit);
    
    %数值计算结果分析和显示
    x = [sol.x];
    y = [sol.y];
    nz = [(sigma_ap/(sigma_ap+sigma_ep)*(y(1,:)+y(2,:))/Ppsat+...
        sigma_as/(sigma_as+sigma_es)*(y(3,:)+y(4,:))/Pssat)./...
        ((y(1,:)+y(2,:))/Ppsat+1+(y(4,:)+y(3,:))/Pssat)];
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
    
%端面抽运的光纤激光器速率方程组
function dy = f(x,y)
    global sigma_ap sigma_ep sigma_as sigma_es gamma_s gamma_p N ...
        alpha_p alpha_s Pssat Ppsat  Ppsp Ppsv
    dy = zeros(4,1); 
    N_2_t=N*(sigma_ap/(sigma_ap+sigma_ep)*(y(1)+y(2))/Ppsat+sigma_as/...
     (sigma_as+sigma_es)*(y(3)+y(4))/Pssat)/((y(1)+y(2))/Ppsat+1+(y(4)+y(3))/Pssat); % 平均上能级粒子数浓度与铒离子平均总掺杂浓度之比
 
    dy(1)=(-gamma_p*(sigma_ap*N-(sigma_ap+sigma_ep)*N21)-alpha_p)*y(1);
    dy(2)=-(-gamma_p*(sigma_ap*N-(sigma_ap+sigma_ep)*N21)-alpha_p)*y(2);
    dy(3)=(gamma_s*((sigma_as+sigma_es)*N21-sigma_as*N)-alpha_s)*y(3);
    dy(4)=-(gamma_s*((sigma_as+sigma_es)*N21-sigma_as*N)-alpha_s)*y(4);

%端面抽运的光纤激光器边界条件
function res = fsbc(y0,yL)
    global R1 R2 Ppl Ppr
    res = [y0(1)-Ppl
       yL(2)-Ppr
       y0(3)-R1*y0(4)
       yL(4)-R2*yL(3)];