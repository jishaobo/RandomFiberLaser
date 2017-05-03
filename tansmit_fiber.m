  %-----------------传输矩阵的功能函数------------------------------?
  function [m_g] =tansmit_fiber(L,kappa,num,n_eff,n_eff0,lambda,lambda_Brag)
  
      delta(num)=2*pi*((n_eff+n_eff0)*1./lambda(num)-n_eff*1./lambda_Brag);
      s(num)=sqrt(kappa.^2-delta(num).^2);

      s11(num)=cosh(s(num)*L)-1i*(delta(num)/s(num))*sinh(s(num)*L);
      s12(num)=-1i*(kappa/s(num))*sinh(s(num)*L);
      s21(num)=1i*(kappa/s(num))*sinh(s(num)*L);
      s22(num)=cosh(s(num)*L)+1i*(delta(num)/s(num))*sinh(s(num)*L);
      m_g=[s11(num) s12(num);s21(num) s22(num)];

  end