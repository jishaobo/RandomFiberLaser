MOIdata = importdata('F:\test.txt');
wv = MOIdata.data(:,1); %����
pw1 = MOIdata.data(:,2); %ͨ��1����
pw2 = MOIdata.data(:,3); %ͨ��2����

figure(1);
subplot(2,1,1);
plot(wv,pw1);
subplot(2,1,2);
plot(wv,pw2);

p_pw1 = max(pw1);
p_idx1 = find(pw1 == p_pw1);
p_wv1 = wv(p_idx1);
disp('ͨ��1��ֵ��Ӧ������');disp(p_wv1);
disp('ͨ��1��ֵ��Ӧ���ʣ�');disp(p_pw1);

p_pw2 = max(pw2);
p_idx2 = find(pw2 == p_pw2);
p_wv2 = wv(p_idx2);
disp('ͨ��2��ֵ��Ӧ������');disp(p_wv2);
disp('ͨ��2��ֵ��Ӧ���ʣ�');disp(p_pw2);