function trueParams = setParamsReduced3state(tres)

% mcmcInfo = struct;
trueParams = struct;
trueParams.discrete_data_flag = true;
trueParams.tres = tres;
trueParams.nStates = 3;

% define transition rate parameters
kon = 0.01;
koff = 0.02;
k_corr_factor = 0.18;
trueParams.R = Q_helper_fun(kon,koff,k_corr_factor);
trueParams.A = expm(trueParams.R*trueParams.tres);                  
trueParams.kon = kon;
trueParams.koff = koff;
trueParams.k_corr_factor = k_corr_factor;
                    
% define emissions vector
r = 0.1;
trueParams.r = r;
r_corr_factor = log(1.3);
trueParams.r_corr_factor = r_corr_factor;

trueParams.v = [0, r, 2*r*exp(r_corr_factor)]'*tres;
trueParams.sigma = 1.5;
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps = 7;
trueParams.alpha_frac = 30/140;
