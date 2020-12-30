% Script to experiment with non-sequential MCMC sampling for parameter
% inference

clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = struct;

%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_mcmc_steps = 50; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.n_chains = 10;

% employ parallel inference?
mcmcInfo.par_chain_flag = true;

% define parameters
mcmcInfo.R = [-.02, .04, 0; .02 -.05 .08; 0 .01 -.08];
mcmcInfo.tres = 20;
mcmcInfo.A = expm(mcmcInfo.R*mcmcInfo.tres);
mcmcInfo.nStates = size(mcmcInfo.A,1);
mcmcInfo.v = [.05,2,4]';
mcmcInfo.seq_length = 60*60/mcmcInfo.tres;

mcmcInfo.nSteps = 7;
[V, D] = eig(mcmcInfo.A);
[~, mi] = max(real(diag(D)));
mcmcInfo.pi0 = V(:,mi)/sum(V(:,mi));

mcmcInfo.sigma = .5;
mcmcInfo.alpha = 0;
mcmcInfo.n_traces = 10;
mcmcInfo.eps = 1e-2;



%%%%%%%%%%%%%%%% Generate helper arrays %%%%%%%%%%%%%%%%
mcmcInfo.coeff_MS2 = ms2_loading_coeff(mcmcInfo.alpha, mcmcInfo.nSteps)';
mcmcInfo.state_options = 1:mcmcInfo.nStates;
mcmcInfo.state_ref = repmat(reshape(mcmcInfo.state_options,1,1,[]),1,mcmcInfo.n_chains);

mcmcInfo.observed_fluo = NaN(mcmcInfo.seq_length,mcmcInfo.n_traces);
mcmcInfo.masterSimStruct = struct;
for n = 1:mcmcInfo.n_traces
    synthetic_data = synthetic_prob(mcmcInfo.seq_length, mcmcInfo.alpha, mcmcInfo.nStates, ...
                          mcmcInfo.nSteps, mcmcInfo.A, mcmcInfo.v', mcmcInfo.sigma, mcmcInfo.pi0);                                     
    
    mcmcInfo.observed_fluo(:,n) = synthetic_data.fluo_MS2;
    % record full simulation info
    fieldNames = fieldnames(synthetic_data);
    for f = 1:length(fieldNames)
        mcmcInfo.masterSimStruct(n).(fieldNames{f}) = synthetic_data.(fieldNames{f});
    end
end

%% %%%%%%%%%%%%%%%%%%%%%% Nonsequential MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize arrays to store inference results
mcmcInfo.logL_vec = NaN(mcmcInfo.n_mcmc_steps,1);
mcmcInfo.A_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_mcmc_steps);
mcmcInfo.v_inf_array = NaN(mcmcInfo.n_mcmc_steps,mcmcInfo.nStates);
mcmcInfo.sigma_inf_array = NaN(mcmcInfo.n_mcmc_steps,1);

mcmcInfo.step = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mcmcInfo = initialize_mcmc_parameters(mcmcInfo);

% mcmcInfo.v_prop_sigma = .1*v2;
% initialize chains
if ~mcmcInfo.par_chain_flag
    mcmcInfo = initialize_chains(mcmcInfo);
else
    mcmcInfo = initialize_chains_par(mcmcInfo);
end

% get predicted fluorescence
mcmcInfo = predict_fluo_full(mcmcInfo,mcmcInfo.v_curr);

wb = waitbar(1,'conducting MCMC inference...');
tic
for step = 2:10%mcmcInfo.n_mcmc_steps
  
    waitbar(step/mcmcInfo.n_mcmc_steps,wb);    
    mcmcInfo.step = step;    
    
    % resample chains   
    tic
    mcmcInfo = resample_chains(mcmcInfo);  
    toc
    % get empirical transition and occupancy counts    
    tic
    mcmcInfo = get_empirical_counts(mcmcInfo);    
    toc
    % use Gibbs sampling to update hyperparameters    
    mcmcInfo = update_hmm_parameters_gibbs(mcmcInfo);
end
toc
disp('done')
delete(wb);