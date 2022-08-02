function jump_weight = calculate_jump_weights2(k1,k2,T)

jump_weight =(k1+(-1).*k2).^(-1).*(k1+(-1).*exp(1).^((-1).*k2.*T).*k1+((-1)+ ...
  exp(1).^((-1).*k1.*T)).*k2);