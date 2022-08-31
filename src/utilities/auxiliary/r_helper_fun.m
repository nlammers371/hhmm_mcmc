function r_vec = r_helper_fun(r,r_corr_factor)

r_corr_factor = exp(r_corr_factor);

r_vec = [0  r  2*r*r_corr_factor];