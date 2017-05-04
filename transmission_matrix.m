  %-----------------��������ܺ���------------------------------
  function [m_g] = transmission_matrix(L,s,n_eff,n_eff0,lambda,num,lambda_Brag)
  
      %L����դ���� 
      %kappa���Ⲩģʽ�����ϵ�� 
      %s�����ƿɼ���
      %n_eff����о��Ч������ 
      %n_eff0�������ʵ������ 
      %lambda������ 
      %lambda_Brag����դ���Ĳ��� 
      kappa(num) = pi*s*n_eff0./lambda(num);
%       delta = 2*pi*((n_eff+n_eff0)*1./lambda-n_eff*1./lambda_Brag);
      delta(num) = 2*pi*n_eff*(1/lambda(num)-1/lambda_Brag);
      omege(num) = sqrt(kappa(num)^2-delta(num)^2); 
      s11(num) = cosh(omege(num)*L)-1i*(delta(num)/omege(num))*sinh(omege(num)*L);
      s12(num) = -1i*(kappa(num)/omege(num))*sinh(omege(num)*L);
      s21(num) = 1i*(kappa(num)/omege(num))*sinh(omege(num)*L);
      s22(num) = cosh(omege(num)*L)+1i*(delta(num)/omege(num))*sinh(omege(num)*L);
      m_g = [s11(num) s12(num);s21(num) s22(num)];

  end