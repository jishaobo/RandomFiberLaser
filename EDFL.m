function  Pout = fiberlaser_twoend
    global R1 R2 Ppl Ppr sigma_ap sigma_ep sigma_as sigma_es qe N Pssat Ppsat L mu h nu_s 
    
    %----------------��������----------------
    lambda_s = 1550 * 1e-9; %
    lambda_p = 976e-9; % ���ֲ�����m��
    tau = 8e-3; % �������ӵ��Է�˥������(����̬����)
    A_c = 3.1416e-10; % ��о�����
    sigma_ap = 1.7e-25; % �������������m^2��
    sigma_ep = 0; % ���ַ��������m^2��
    [sigma_as,sigma_es] = GetErSpectrum(lambda_s);
    N = 6e25; % ��λ���������������m^-3��
    qe = 0.85; % ��������Ч��
    R1 =.99;
    R2 =.035;
    L = 40;
    mu = 2; % ƫ��ģ��Ŀ
    
    %----------------���������м���̲���----------------
    c = 3e8; % ����й���
    h = 6.626e-34; % ���ʿ˳���
    nu_s = c/lambda_s; % 
    nu_p = c/lambda_p; % ���ֹ�Ƶ��   
    Pth = h * nu_s * A_c * N / (tau * pi); % ��ֵ����
    Pssat = h * nu_s * A_c * N /( qe * (sigma_es+sigma_as) * tau); % ���͹���
    Ppsat = h * nu_p * A_c * N /( qe * (sigma_ep+sigma_ap) * tau);
    
    %���˹⹦������
    Ppl = 100;
    Ppr = 0;
    
    %----------------��ֵ������ֵ���----------------
    OPTION = bvpset('Stats','ON');
    solinit = bvpinit(linspace(0,L,100),[Ppl Ppr 30 Ppr]);
    sol = bvp4c(@f,@fsbc,solinit,OPTION);
    
    %----------------��ֵ��������������ʾ----------------
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
   
    
%----------------��Ŀ�겨��������˼��ⷢ��/�������----------------    
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
    abint = interp1(wl,ab,wlint,'PCHIP');%�ֶ�����Hermite��ֵ
    emint = interp1(wl,em,wlint,'PCHIP');
    
    
%----------------��Ŀ�겨��������˼��ⷢ��/�������----------------
function dy = f(x,y)
    global sigma_ap sigma_ep sigma_as sigma_es qe N Pssat Ppsat mu h nu_s 
    delta_nu_s = ;
    dy = zeros(3,1); 
    Nz=N*(sigma_ap/(sigma_ap+sigma_ep)*(y(1)+y(2))/Ppsat+sigma_as/...
     (sigma_as+sigma_es)*(y(3)+y(4))/Pssat)/((y(1)+y(2))/Ppsat+1+(y(4)+y(3))/Pssat); % ƽ�����ܼ�������Ũ����������ƽ���ܲ���Ũ��֮��
 
    dy(1) = (qe*(sigma_ap+sigma_ep)*Nz - qe*sigma_ap) * y(1);
    dy(2) = (qe*(sigma_ap+sigma_ep)*Nz - qe*sigma_ap) * y(2) + mu * qe*sigma_ep*Nz*h*nu_s*delta_nu_s;
    dy(3) = qe*sigma_ap*Nz * y(3);

    
%----------------�߽�����----------------
function res = fsbc(y0,yL)
    global R1 R2 Ppl Ppr
    res = [y0(1)-Ppl
       yL(2)-Ppr
       y0(3)-R1*y0(4)
       yL(4)-R2*yL(3)];