function loss = calculateModelLoss(param_vec)

global mcmcInfo

n_chains = mcmcInfo.n_chains;

mcmcInfo.k_curr = param_vec(1:3);
kc = mcmcInfo.k_curr;
kc(end) = exp(kc(end));
mcmcInfo.r_curr = param_vec(4:5);
rc = mcmcInfo.r_curr;
rc(end) = exp(rc(end));
mcmcInfo.sigma_curr = repmat(param_vec(end),n_chains,1);
% sg = mcmcInfo.sigma_curr;

% generate full variable arrays
Q_init = Q_helper_fun(kc(1),kc(2),kc(3));
mcmcInfo.A_curr = repmat(expm(Q_init*mcmcInfo.tres),1,1,n_chains);
r_vec = [0 rc(1) 2*rc(1)*rc(2)];

% calculate pi0 
[V, D] = eig(mcmcInfo.A_curr(:,:,1));
[~, mi] = max(real(diag(D)));
mcmcInfo.pi0_curr = repmat((V(:,mi)/sum(V(:,mi)))',n_chains,1);    

mcmcInfo.v_curr = repmat(r_vec*mcmcInfo.tres,n_chains,1);
mcmcInfo = resample_chains_v5(mcmcInfo);    % "Expectation Step"                 

% get empirical transition and occupancy counts    
% mcmcInfo = get_empirical_counts_v3(mcmcInfo); 

% calculate updated logL
mcmcInfo = calculateLogLikelihood(mcmcInfo);

loss = max(mcmcInfo.logL_vec(mcmcInfo.step,:),[],2);

mcmcInfo.step = mcmcInfo.step + 1;