MOIdata = importdata('F:\test.txt');
wv = MOIdata.data(:,1); %波长
pw1 = MOIdata.data(:,2); %通道1功率
pw2 = MOIdata.data(:,3); %通道2功率

figure(1);
subplot(2,1,1);
plot(wv,pw1);
subplot(2,1,2);
plot(wv,pw2);

p_pw1 = max(pw1);
p_idx1 = find(pw1 == p_pw1);
p_wv1 = wv(p_idx1);
disp('通道1峰值对应波长：');disp(p_wv1);
disp('通道1峰值对应功率：');disp(p_pw1);

p_pw2 = max(pw2);
p_idx2 = find(pw2 == p_pw2);
p_wv2 = wv(p_idx2);
disp('通道2峰值对应波长：');disp(p_wv2);
disp('通道2峰值对应功率：');disp(p_pw2);