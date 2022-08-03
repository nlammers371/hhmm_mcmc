function jump_weight = calculate_jump_weights4Rep2(k1,k2,k3,T)

% weight calculation for path of type 1->2->3->2->3
% k1 is the first transition rate
% k2 is the repeated rate
% k3 is the third rate in the sequence

jump_weight = 1+exp(1).^((-1).*k1.*T).*(k1+(-1).*k2).^(-2).*k2.^2.*(k1+(-1).*k3) ...
  .^(-1).*k3+exp(1).^((-1).*k3.*T).*k1.*k2.^2.*(k2+(-1).*k3).^(-2).* ...
  ((-1).*k1+k3).^(-1)+exp(1).^((-1).*k2.*T).*k1.*(k1+(-1).*k2).^(-2) ...
  .*(k2+(-1).*k3).^(-2).*k3.*(2.*k1.*k2+(-3).*k2.^2+(-1).*k1.*k3+2.* ...
  k2.*k3+(k1+(-1).*k2).*k2.*(k2+(-1).*k3).*T);

