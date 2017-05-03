clear all;
clc;
%-------------��������-------------------------------         
n_eff = 1.4452; %n_eff----------------��դ��Ч������
L = 0.1*1e-2; %L----------------��դ����(cm)
lambda_Brag = 1550*1e-9; %lambda_Brag----------��դ���Ĳ���
n_eff0 = 5*1e-6; %n_eff0-----------�����ʵ������
s = 1;%s--------------------�����ʵ��Ƶ����ƿɼ���
k = 2*pi*n_eff/lambda_Brag; %-----------���ϵ��
%kappa_L----------------�Ⲩģʽ�����ϵ��
%delta----------------ģ��ʧг��
%-------------��դ����-------------------------------
disp('����Ĳ�������>>>>>>>>>>>>>>');
disp('��դ��Ч������Ϊ��');disp(n_eff);
disp('��դ�������Ϊ��');disp(n_eff0);
disp('�������դ�����Ĳ���Ϊ');disp(lambda_Brag);
disp('�������դ�ĳ���Ϊ��');disp(L);
lambda = 1e-9*linspace(1548,1552,1000);
F = [1 0;0 1]; %��ʼ����
spacing = (rand(1,1000)+0.5)*L; %������ 0.5��դ��-1.5��դ��
refspacing = ones(1,1000);
%       
%       m_d = [exp(1i*k*spacing(num)) 0;0 exp(-1i*k*spacing(num))]; %����������
% 
%       f = m_g*m_d;
%-------------������դ---------------
for num=1:1000
    kappa=pi*s*n_eff0./lambda(num);
    m_g=tansmit_fiber(L,kappa,num,n_eff,n_eff0,lambda,lambda_Brag);
    m_g=m_g*F;    
    r(num)=m_g(2,1)/m_g(1,1);%������դ���������
    R(num)=(abs(-r(num)))^2;%������դ��ǿ������
end
figure(1);
subplot(2,1,1);
logR = 10*log10(R);
plot(lambda*1e9,logR);
xlim([1548,1552]);
%------------��դ����----------------
m_total = [1 0;0 1];
for idx = 1:1000
    m_d = [exp(1i*k*spacing(idx)) 0;0 exp(-1i*k*spacing(idx))]; %����������
    m_temp = m_g*m_d;
    m_total = m_temp*m_total;
end
finalr = m_total(2,1)/m_total(1,1);%��դ�������������
finalR = (abs(-finalr(num)))^2;%��դ���й�ǿ������
subplot(2,1,2);
logfinalR = 10*log10(finalR);
plot(lambda*1e9,logfinalR);
xlim([1548,1552]);
%-------------------------3dB����-------------------------       
maxData = max(logR);
threshold = maxData- 3;
bandIdx = find(logR>threshold);
startIdx = bandIdx(1);
stopIdx = bandIdx(end);
disp('3dB����'),disp((lambda(stopIdx)-lambda(startIdx))*1e9);

