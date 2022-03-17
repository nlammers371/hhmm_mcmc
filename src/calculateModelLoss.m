function loss = calculateModelLoss(param_vec)

global mcmcInfo

mcmcInfo.v_curr = param_vec(1:3);
mcmcInfo.A_curr = reshape(param_vec(4:12),3,3);
mcmcInfo.A_curr = mcmcInfo.A_curr ./ sum(mcmcInfo.A_curr);
mcmcInfo.sigma_curr = param_vec(end);

mcmcInfo = resample_chains_v3(mcmcInfo);    % "Expectation Step"                 

% get empirical transition and occupancy counts    
mcmcInfo = get_empirical_counts_v3(mcmcInfo);

% use Gibbs sampling to update hyperparameters      
% mcmcInfo = update_hmm_parameters_v3(mcmcInfo);   

% calculate updated logL
mcmcInfo = calculateLogLikelihood(mcmcInfo);

loss = -mcmcInfo.logL_vec(mcmcInfo.step,1);

mcmcInfo.step = mcmcInfo.step + 1;