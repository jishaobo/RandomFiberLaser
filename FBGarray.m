clear all;
clc;
%-------------变量定义-------------------------------  
pointnum = 1000; %设定光谱点数
secnum = 100; %阵列中光栅个数
lambda = 1e-9*linspace(1545,1555,pointnum);
n_eff = 1.4452; %n_eff---------------纤芯有效折射率
L = 0.1*1e-2; %L-----------------------光栅长度(cm)
lambda_Brag = 1550*1e-9; %lambda_Brag--光栅中心波长
n_eff0 = 5*1e-6; %n_eff0-------------折射率调制深度
s = 1;%s----------------------折射率调制的条纹可见度
k = 2*pi*n_eff./lambda; %k-----------------耦合系数
disp('输入的参数如下>>>>>>>>>>>>>>');
disp('光栅有效折射率为：');disp(n_eff);
disp('光栅调制深度为：');disp(n_eff0);
disp('布拉格光栅的中心波长为');disp(lambda_Brag);
disp('布拉格光栅的长度为：');disp(L);
%==============光栅仿真==============
spacing = (rand(1,secnum)+0.5)*L; %随机间隔 0.5倍栅长-1.5倍栅长
refspacing = ones(1,secnum)*L;
%--------------单个光栅---------------
for a = 1:pointnum
    m_g = transmission_matrix(L,s,n_eff,n_eff0,lambda,a,lambda_Brag);  
    r(a) = m_g(2,1)/m_g(1,1);%单个光栅振幅反射率
    R(a) = (abs(r(a)))^2;%单个光栅光强反射率
end
figure(1);
subplot(1,3,1);
logR = 10*log10(R);
plot(lambda*1e9,logR);
xlim([1545,1555]);
title('单个光栅')
%--------------3dB带宽--------------       
maxData = max(logR);
threshold = maxData- 3;
bandIdx = find(logR>threshold);
startIdx = bandIdx(1);
stopIdx = bandIdx(end);
disp('3dB带宽：'),disp((lambda(stopIdx)-lambda(startIdx))*1e9);
%------------随机分布光栅阵列----------------
for a = 1:pointnum
    for b = 1:secnum
        m_g = transmission_matrix(L,s,n_eff,n_eff0,lambda,a,lambda_Brag);  
        m_d = [exp(-1i*k(a)*spacing(b)) 0;0 exp(1i*k(a)*spacing(b))]; %间隔传输矩阵
        if b == 1
            m_total = m_d*m_g;
        else
            m_total = m_d*m_g*m_total;
        end
        finalr(a) = m_total(2,1)/m_total(1,1);%单个光栅振幅反射率
        finalR(a) = (abs(finalr(a)))^2;%单个光栅光强反射率
    end
end
subplot(1,3,2);
logfinalR = 10*log10(finalR);
plot(lambda*1e9,logfinalR);
xlim([1545,1555]);
title('随机分布光栅阵列')
% %------------均匀分布光栅阵列----------------
% for a = 1:pointnum
%     for b = 1:secnum
%         m_g = transmission_matrix(L,s,n_eff,n_eff0,lambda,a,lambda_Brag);  
%         m_d = [exp(-1i*k(a)*refspacing(b)) 0;0 exp(1i*k(a)*refspacing(b))]; %间隔传输矩阵
%         if b == 1
%             m_total = m_d*m_g;
%         else
%             m_total = m_d*m_g*m_total;
%         end
%         refr(a) = m_total(2,1)/m_total(1,1);%单个光栅振幅反射率
%         refR(a) = (abs(refr(a)))^2;%单个光栅光强反射率
%     end
% end
% subplot(1,3,3);
% logrefR = 10*log10(refR);
% plot(lambda*1e9,logrefR);
% xlim([1545,1555]);
% title('均匀分布光栅阵列')

%------------调整数据结构----------------
lambda = (lambda*10^9)';
logR = logR';
logfinalR = logfinalR';
% logrefR = logrefR';

