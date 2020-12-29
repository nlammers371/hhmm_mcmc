% Script to experiment with non-sequential MCMC sampling for parameter
% inference

clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = struct;

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

%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_mcmc_steps = 100; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.n_chains = 100;

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

%%% Nonsequential MCMC

% initialize arrays to store inference results
mcmcInfo.logL_vec = NaN(mcmcInfo.n_mcmc_steps,1);
mcmcInfo.A_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_mcmc_steps);
mcmcInfo.v_inf_array = NaN(mcmcInfo.n_mcmc_steps,mcmcInfo.nStates);
mcmcInfo.sigma_inf_array = NaN(mcmcInfo.n_mcmc_steps,1);

% initialize arrays to store inference diagnostics
mcmcInfo.v_acceptance_array = NaN(mcmcInfo.n_mcmc_steps,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A prior--assume strongly diagonal PDF given short timescale
% take A columns to follow multinomial Dirichlet distribution
mcmcInfo.A_alpha = ones(mcmcInfo.nStates);%*n_particles*n_traces;
mcmcInfo.A_alpha(eye(mcmcInfo.nStates)==1) = mcmcInfo.A_alpha(eye(mcmcInfo.nStates)==1)*10; % distribution hyper params
mcmcInfo.A_curr = mcmcInfo.A;%sample_A_dirichlet(mcmcInfo.A_alpha, zeros(mcmcInfo.nStates));

% calculate pi0 
[V, D] = eig(mcmcInfo.A_curr);
[~, mi] = max(real(diag(D)));
mcmcInfo.pi0_curr = V(:,mi)/sum(V(:,mi));

% initialize sigma as inverse gamma (see: http://ljwolf.org/teaching/gibbs.html)
fluo_vec = mcmcInfo.observed_fluo(:);
mcmcInfo.sigma_curr = mcmcInfo.sigma;%0.1*mean(fluo_vec);%sqrt(1./gamrnd(100*mcmcInfo.seq_length*mcmcInfo.n_traces/2,1./(fluo_vec'*fluo_vec)));

% initialize v
v2 = prctile(fluo_vec,95) / 7;%mean(fluo_vec)/sum(mcmcInfo.coeff_MS2)/(mcmcInfo.pi0_curr(2)+2*mcmcInfo.pi0_curr(3));
mcmcInfo.v_curr = mcmcInfo.v;%[0 ; v2 ; 2*v2] + rand(3,1)*.2;
% mcmcInfo.v_prop_sigma = .1*v2;

% initialize chains
mcmcInfo = initialize_chains(mcmcInfo);

% get predicted fluorescence
mcmcInfo = predict_fluo_full(mcmcInfo);

wb = waitbar(0,'conducting MCMC inference...');
for step = 1:50%mcmcInfo.n_mcmc_steps %mcmcInfo.n_mcmc_steps    
    waitbar(step/mcmcInfo.n_mcmc_steps,wb);
    
    mcmcInfo.step = step;
    
    % resample chains
    mcmcInfo = resample_chains(mcmcInfo);          
    
    % get empirical transition and occupancy counts
    mcmcInfo = get_empirical_counts(mcmcInfo);
    
    % use Gibbs sampling to update hyperparameters
    mcmcInfo = update_hmm_parameters_v1(mcmcInfo);
    
end
disp('done')
delete(wb);