clear all;
clc;
%-------------��������-------------------------------  
pointnum = 1000; %�趨���׵���
secnum = 100; %�����й�դ����
lambda = 1e-9*linspace(1545,1555,pointnum);
n_eff = 1.4452; %n_eff---------------��о��Ч������
L = 0.1*1e-2; %L-----------------------��դ����(cm)
lambda_Brag = 1550*1e-9; %lambda_Brag--��դ���Ĳ���
n_eff0 = 5*1e-6; %n_eff0-------------�����ʵ������
s = 1;%s----------------------�����ʵ��Ƶ����ƿɼ���
k = 2*pi*n_eff./lambda; %k-----------------���ϵ��
disp('����Ĳ�������>>>>>>>>>>>>>>');
disp('��դ��Ч������Ϊ��');disp(n_eff);
disp('��դ�������Ϊ��');disp(n_eff0);
disp('�������դ�����Ĳ���Ϊ');disp(lambda_Brag);
disp('�������դ�ĳ���Ϊ��');disp(L);
%==============��դ����==============
spacing = (rand(1,secnum)+0.5)*L; %������ 0.5��դ��-1.5��դ��
refspacing = ones(1,secnum)*L;
%--------------������դ---------------
for a = 1:pointnum
    m_g = transmission_matrix(L,s,n_eff,n_eff0,lambda,a,lambda_Brag);  
    r(a) = m_g(2,1)/m_g(1,1);%������դ���������
    R(a) = (abs(r(a)))^2;%������դ��ǿ������
end
figure(1);
subplot(1,3,1);
logR = 10*log10(R);
plot(lambda*1e9,logR);
xlim([1545,1555]);
title('������դ')
%--------------3dB����--------------       
maxData = max(logR);
threshold = maxData- 3;
bandIdx = find(logR>threshold);
startIdx = bandIdx(1);
stopIdx = bandIdx(end);
disp('3dB����'),disp((lambda(stopIdx)-lambda(startIdx))*1e9);
%------------����ֲ���դ����----------------
for a = 1:pointnum
    for b = 1:secnum
        m_g = transmission_matrix(L,s,n_eff,n_eff0,lambda,a,lambda_Brag);  
        m_d = [exp(-1i*k(a)*spacing(b)) 0;0 exp(1i*k(a)*spacing(b))]; %����������
        if b == 1
            m_total = m_d*m_g;
        else
            m_total = m_d*m_g*m_total;
        end
        finalr(a) = m_total(2,1)/m_total(1,1);%������դ���������
        finalR(a) = (abs(finalr(a)))^2;%������դ��ǿ������
    end
end
subplot(1,3,2);
logfinalR = 10*log10(finalR);
plot(lambda*1e9,logfinalR);
xlim([1545,1555]);
title('����ֲ���դ����')
% %------------���ȷֲ���դ����----------------
% for a = 1:pointnum
%     for b = 1:secnum
%         m_g = transmission_matrix(L,s,n_eff,n_eff0,lambda,a,lambda_Brag);  
%         m_d = [exp(-1i*k(a)*refspacing(b)) 0;0 exp(1i*k(a)*refspacing(b))]; %����������
%         if b == 1
%             m_total = m_d*m_g;
%         else
%             m_total = m_d*m_g*m_total;
%         end
%         refr(a) = m_total(2,1)/m_total(1,1);%������դ���������
%         refR(a) = (abs(refr(a)))^2;%������դ��ǿ������
%     end
% end
% subplot(1,3,3);
% logrefR = 10*log10(refR);
% plot(lambda*1e9,logrefR);
% xlim([1545,1555]);
% title('���ȷֲ���դ����')

%------------�������ݽṹ----------------
lambda = (lambda*10^9)';
logR = logR';
logfinalR = logfinalR';
% logrefR = logrefR';

