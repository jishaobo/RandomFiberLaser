clear;
clc;

%----------------变量设定----------------
n = 1.4682; %纤芯有效折射率
L_array = [200e3, 50e3, 10e3]; %延时光纤长度
Delta_nu = 1e3; %激光器真实线宽
alpha = 1; %耦合器分光比
Omega = 200e6; %AOM移频量
c = 3e8; %真空中                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              光速
I_0 = 1; %入射光光强
omega = 199.5e6:0.1:200.5e6; %功率谱角频率
tau_c = 0.318/Delta_nu; %激光器相干时间

figure(1);
cm = colormap(jet(length(L_array))); 
for a = 1:length(L_array)
    L = L_array(a);
    tau_d = L*n/c; %光纤延时线对应的延迟时间
    
    PSD = 0.5*alpha^2*I_0^2* (  2/tau_c./( (2/tau_c)^2 + (omega-Omega).^2 ) .* ( 1-exp(-2*tau_d/tau_c).*( cos((omega-Omega)*tau_d) +...
    +2/tau_c./(omega-Omega).*sin((omega-Omega)*tau_d) ) ) + pi*exp(-2*tau_d/tau_c).*dirac(omega-Omega) );

    log_PSD = 10*log10(PSD);

    plot(omega,log_PSD,'color',cm(a,:));
    leg_r{a} = ['延时光纤长度 ',num2str(L)];
    hold on;
end
legend(leg_r);

