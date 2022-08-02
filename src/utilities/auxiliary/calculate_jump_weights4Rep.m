function jump_weight = calculate_jump_weights4Rep(k1,k2,T)

jump_weight = (k1+(-1).*k2).^(-3).*((k1+(-1).*k2).^3+(-1).*exp(1).^((-1).*k2.*T) ...
  .*k1.^2.*(k1+k1.*k2.*T+(-1).*k2.*(3+k2.*T))+exp(1).^((-1).*k1.*T) ...
  .*k2.^2.*(k2+k1.*((-3)+(-1).*k1.*T+k2.*T)));

  

  
