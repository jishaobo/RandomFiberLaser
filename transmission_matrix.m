  %-----------------传输矩阵功能函数------------------------------
  function [m_g] = transmission_matrix(L,s,n_eff,n_eff0,lambda,lambda_Brag)
  
      %L：光栅长度 
      %kappa：光波模式互耦合系数 
      %s：条纹可见度
      %n_eff：纤芯有效折射率 
      %n_eff0：折射率调制深度 
      %lambda：波长 
      %lambda_Brag：光栅中心波长 
      kappa = pi*s*n_eff0./lambda;
%       delta = 2*pi*((n_eff+n_eff0)*1./lambda-n_eff*1./lambda_Brag);
      delta = 2*pi*n_eff*(1./lambda-1/lambda_Brag);
      omege = sqrt(kappa.^2-delta.^2); 
      s11 = cosh(omege*L)-1i*(delta./omege).*sinh(omege*L);
      s12 = -1i*(kappa./omege).*sinh(omege*L);
      s21 = 1i*(kappa./omege).*sinh(omege*L);
      s22 = cosh(omege*L)+1i*(delta./omege).*sinh(omege*L);
      m_g=[s11 s12;s21 s22];

  end