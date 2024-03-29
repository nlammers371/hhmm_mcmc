function jump_weight = calculate_jump_weights3(k1,k2,k3,T)

jump_weight = 1+(k1+(-1).*k2).^(-1).*(k1+(-1).*k3).^(-1).*(k2+(-1).*k3).^(-1).*( ...
  exp(1).^((-1).*k3.*T).*k1.*k2.*((-1).*k1+k2)+exp(1).^((-1).*k2.*T) ...
  .*k1.*(k1+(-1).*k3).*k3+exp(1).^((-1).*k1.*T).*k2.*k3.*((-1).*k2+ ...
  k3));