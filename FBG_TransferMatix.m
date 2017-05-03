clear all;
clc;
%-------------变量定义-------------------------------         
n_eff = 1.4452; %n_eff----------------光栅有效折射率
L = 0.1*1e-2; %L----------------光栅长度(cm)
lambda_Brag = 1550*1e-9; %lambda_Brag----------光栅中心波长
n_eff0 = 5*1e-6; %n_eff0-----------折射率调制深度
s = 1;%s--------------------折射率调制的条纹可见度
k = 2*pi*n_eff/lambda_Brag; %-----------耦合系数
%kappa_L----------------光波模式互耦合系数
%delta----------------模间失谐量
%-------------光栅仿真-------------------------------
disp('输入的参数如下>>>>>>>>>>>>>>');
disp('光栅有效折射率为：');disp(n_eff);
disp('光栅调制深度为：');disp(n_eff0);
disp('布拉格光栅的中心波长为');disp(lambda_Brag);
disp('布拉格光栅的长度为：');disp(L);
lambda = 1e-9*linspace(1548,1552,1000);
F = [1 0;0 1]; %初始条件
spacing = (rand(1,1000)+0.5)*L; %随机间隔 0.5倍栅长-1.5倍栅长
refspacing = ones(1,1000);
%       
%       m_d = [exp(1i*k*spacing(num)) 0;0 exp(-1i*k*spacing(num))]; %间隔传输矩阵
% 
%       f = m_g*m_d;
%-------------单个光栅---------------
for num=1:1000
    kappa=pi*s*n_eff0./lambda(num);
    m_g=tansmit_fiber(L,kappa,num,n_eff,n_eff0,lambda,lambda_Brag);
    m_g=m_g*F;    
    r(num)=m_g(2,1)/m_g(1,1);%单个光栅振幅反射率
    R(num)=(abs(-r(num)))^2;%单个光栅光强反射率
end
figure(1);
subplot(2,1,1);
logR = 10*log10(R);
plot(lambda*1e9,logR);
xlim([1548,1552]);
%------------光栅阵列----------------
m_total = [1 0;0 1];
for idx = 1:1000
    m_d = [exp(1i*k*spacing(idx)) 0;0 exp(-1i*k*spacing(idx))]; %间隔传输矩阵
    m_temp = m_g*m_d;
    m_total = m_temp*m_total;
end
finalr = m_total(2,1)/m_total(1,1);%光栅阵列振幅反射率
finalR = (abs(-finalr(num)))^2;%光栅阵列光强反射率
subplot(2,1,2);
logfinalR = 10*log10(finalR);
plot(lambda*1e9,logfinalR);
xlim([1548,1552]);
%-------------------------3dB带宽-------------------------       
maxData = max(logR);
threshold = maxData- 3;
bandIdx = find(logR>threshold);
startIdx = bandIdx(1);
stopIdx = bandIdx(end);
disp('3dB带宽：'),disp((lambda(stopIdx)-lambda(startIdx))*1e9);

