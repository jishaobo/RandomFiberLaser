%--------------------CSV格式读入---------------------------
files=dir('F:\FBG\实验\20170522_tunable\可调滤波器\一组\*.*.csv');
pwlist = [];
maxwvlist = [];
maxpwlist = [];
disp(length(files));
for i = 1:length(files)
    f=fullfile('F:','FBG','实验','20170522_tunable','可调滤波器','一组',files(i).name);
    disp(f);
    ex = csvread(f,29,0);
    wv = ex(:,1);
    pw = ex(:,2);   
    maxpw = max(pw);
    maxidx = find(pw==maxpw);
    maxwv = wv(maxidx(1));  %避免出现多个值导致后面矩阵大小不一致
    pwlist = [pwlist,pw];
    maxwvlist = [maxwvlist,maxwv];
    maxpwlist = [maxpwlist,maxpw];
end
maxwvlist = maxwvlist';
deltawv = max(maxwvlist) - min(maxwvlist);
disp('波长波动为：');
disp(deltawv);
maxpwlist = maxpwlist';
deltapw = max(maxpwlist) - min(maxpwlist);
disp('功率波动为：');
disp(deltapw);
