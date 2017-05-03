  %-----------------��������ܺ���------------------------------
  function [m_g] = transmission_matrix(L,s,n_eff,n_eff0,lambda,lambda_Brag)
  
      %L����դ���� 
      %kappa���Ⲩģʽ�����ϵ�� 
      %s�����ƿɼ���
      %n_eff����о��Ч������ 
      %n_eff0�������ʵ������ 
      %lambda������ 
      %lambda_Brag����դ���Ĳ��� 
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