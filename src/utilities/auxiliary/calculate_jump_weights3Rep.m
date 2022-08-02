function jump_weight = calculate_jump_weights3Rep(k1,k2,T)

jump_weight = exp(1).^((-1).*(k1+k2).*T).*(k1+(-1).*k2).^(-2).*((-1).*exp(1).^( ...
  k1.*T).*k1.^2+exp(1).^((k1+k2).*T).*(k1+(-1).*k2).^2+(-1).*exp(1) ...
  .^(k2.*T).*k2.*(k2+k1.*((-2)+(-1).*k1.*T+k2.*T)));


  

  
