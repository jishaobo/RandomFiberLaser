%--------------------CSV��ʽ����---------------------------
files=dir('F:\FBG\ʵ��\20170522_tunable\�ɵ��˲���\һ��\*.*.csv');
pwlist = [];
maxwvlist = [];
maxpwlist = [];
disp(length(files));
for i = 1:length(files)
    f=fullfile('F:','FBG','ʵ��','20170522_tunable','�ɵ��˲���','һ��',files(i).name);
    disp(f);
    ex = csvread(f,29,0);
    wv = ex(:,1);
    pw = ex(:,2);   
    maxpw = max(pw);
    maxidx = find(pw==maxpw);
    maxwv = wv(maxidx(1));  %������ֶ��ֵ���º�������С��һ��
    pwlist = [pwlist,pw];
    maxwvlist = [maxwvlist,maxwv];
    maxpwlist = [maxpwlist,maxpw];
end
maxwvlist = maxwvlist';
deltawv = max(maxwvlist) - min(maxwvlist);
disp('��������Ϊ��');
disp(deltawv);
maxpwlist = maxpwlist';
deltapw = max(maxpwlist) - min(maxpwlist);
disp('���ʲ���Ϊ��');
disp(deltapw);
