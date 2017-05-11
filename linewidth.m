clear;
clc;

%----------------�����趨----------------
n = 1.4682; %��о��Ч������
L_array = [200e3, 50e3, 10e3]; %��ʱ���˳���
Delta_nu = 1e3; %��������ʵ�߿�
alpha = 1; %������ֹ��
Omega = 200e6; %AOM��Ƶ��
c = 3e8; %�����                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ����
I_0 = 1; %������ǿ
omega = 199.5e6:0.1:200.5e6; %�����׽�Ƶ��
tau_c = 0.318/Delta_nu; %���������ʱ��

figure(1);
cm = colormap(jet(length(L_array))); 
for a = 1:length(L_array)
    L = L_array(a);
    tau_d = L*n/c; %������ʱ�߶�Ӧ���ӳ�ʱ��
    
    PSD = 0.5*alpha^2*I_0^2* (  2/tau_c./( (2/tau_c)^2 + (omega-Omega).^2 ) .* ( 1-exp(-2*tau_d/tau_c).*( cos((omega-Omega)*tau_d) +...
    +2/tau_c./(omega-Omega).*sin((omega-Omega)*tau_d) ) ) + pi*exp(-2*tau_d/tau_c).*dirac(omega-Omega) );

    log_PSD = 10*log10(PSD);

    plot(omega,log_PSD,'color',cm(a,:));
    leg_r{a} = ['��ʱ���˳��� ',num2str(L)];
    hold on;
end
legend(leg_r);

