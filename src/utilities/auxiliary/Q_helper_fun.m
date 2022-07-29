function Q_init = Q_helper_fun(kon,koff,k_corr_factor)

k_corr_factor = exp(k_corr_factor);

Q_init = [-2*kon    koff*k_corr_factor         0;
           2*kon   -(koff+kon)*k_corr_factor   2*koff;
               0    kon*k_corr_factor         -2*koff]; 