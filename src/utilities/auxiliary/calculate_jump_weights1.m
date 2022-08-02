function jump_weight = calculate_jump_weights1(k1,T)

jump_weight = 1+(-1).*exp(1).^((-1).*k1.*T);
