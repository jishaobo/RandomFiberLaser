function G = EDFAgain(wl,L,Pp)

    % wl 输入波长
    % L EDF长度
    % Pp 泵浦功率
    
    %-----------------参数设置-----------------
    wlp = 976e-9; % 泵浦波长（m）
    loss = 2; % 系统总损耗（2dB）

    d = 6e-6; % 掺铒光纤模场直径（m）
    Gamma = 0.722; 
    A = pi*((d/2)^2)/Gamma;  % 掺铒光纤有效面积（m^2）
    h = 6.626068e-034; % 普朗克常量
    c = 3e8; % 真空中光速
    nup = c/wlp; % 泵浦光频率
    doping = 500; % 分子分数
    
    % -----------------关于掺铒光纤-----------------
    sap = 1.7e-25; % 泵浦吸收面积（m^2）
    sep = 0; % 泵浦发射面积（m^2）
    [sel,sal] = GetErSpectrum(wl);
    N = 6e25; % 单位体积掺杂离子数（m^-3）
    tau = 8e-3; % 激发离子的自发衰变速率
    qe = 0.85; % 泵浦量子效率
    fiberStr = 'Er-Doped';

    % -----------------饱和功率-----------------
    PpSat = (h*nup*A)/((sap+sep)*tau*qe); % 泵浦饱和功率

    % -----------------计算单次增益-----------------
    % 1. 计算吸收功率
    [Ppl,err] = SolveForPpz(Pp,PpSat,-N*sap*L);
    Pa = Pp-Ppl;
    % 2. 计算增益（dB）
    for ww = 1:numel(wl)
         (ww,Pp,L) = 4.34*((qe*(sal(ww)+sel(ww))*tau*Pa)./(A*h*nup) - Gamma*N*sal(ww)*L) - loss; % 各波长处的增益
    end
    % 3. 定义图例
    legStr{Pp} = sprintf('P_p = %d mW',round(1000*Pp));
       
    % -----------------画图-----------------
    cmap = hsv(length(L));
    figure('windowstyle','docked');
    for L = 1:length(L);
        plot(1000*Pp,G(1,:,L),'linewidth',2,'color',cmap(L,:)); hold on;
        legStr{L} = sprintf('%.2f m',L(L));
    end
    hT = title(sprintf('Single-Pass Gain at %d nm\nUsing %s Fiber (%d ppm)',round(wl*1e9),fiberStr,10*round(doping/10))); set(hT,'fontsize',16);
    set(gca,'fontsize',14,'linewidth',2); grid on;
    xlabel(sprintf('Pump power (mW) at %d nm',round(wlp*1e9))); ylabel('Single-Pass Gain (dB)');
    legend(legStr,'Location','Best');        

end



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
    abint = interp1(wl,ab,wlint,'cubic');
    emint = interp1(wl,em,wlint,'cubic');

end



function [solution,residual] = SolveForPpz(Pp0,Ps,k)
    % based on equation: log(Pp(z) / Pp(0)) + (Pp(z) - Pp(0))/Ps = -N*sig_ap*a
    % where k is the right hand side (should be negative)

    if nargin<1
        Pp0 = 0.01;
        Ps = 0.002; % about right for 915 nm
        N = 1.5e25; % 550 ppm doping
        sap = 8e-25; % about right for 915 nm
        z = 1;
        k = -N*sap*z;
    end

    % initial guess is that between 0.01 and 0.99 of the power is absorbed
    Ppz = Pp0*exp(k);
    res = log(Ppz/Pp0) + (Ppz-Pp0)/Ps - k;
    % solution = Ppz;
    % residual = res;

    % solving numerically, we can incorporate the saturation term as well
    absDb = [0:0.001:30];
    txLin = 10.^(-absDb./10);
    Ppz = txLin*Pp0; 

    % calculate residual for each Ppz guess
    res = log(Ppz./Pp0) + (Ppz-Pp0)./Ps - k;

    % output value with lowest residual
    [minres,iimin] = min(abs(res));
    solution = Ppz(iimin);
    residual = minres;

end