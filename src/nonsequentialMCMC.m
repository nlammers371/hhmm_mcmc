% Script to experiment with non-sequential MCMC sampling for parameter
% inference

clear 
close all

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = struct;

% define parameters
mcmcInfo.R = [-.02, .04, 0; .02 -.05 .08; 0 .01 -.08];
mcmcInfo.tres = 20;
mcmcInfo.A = expm(mcmcInfo.R*mcmcInfo.tres);
mcmcInfo.nStates = size(mcmcInfo.A,1);
mcmcInfo.v = [.05,2,4];
mcmcInfo.seq_length = 60*60/mcmcInfo.tres;
mcmcInfo.nSteps = 7;
[V, D] = eig(mcmcInfo.A);
[~, mi] = max(real(diag(D)));
mcmcInfo.pi0 = V(:,mi)/sum(V(:,mi));

mcmcInfo.noise = 1;
mcmcInfo.alpha = 1.4;
mcmcInfo.n_traces = 10;
mcmcInfo.eps = 1e-2;

mcmcInfo.coeff_MS2 = ms2_loading_coeff(mcmcInfo.alpha, mcmcInfo.nSteps)';

mcmcInfo.masterSimStruct = struct;
for i = 1:mcmcInfo.n_traces
    synthetic_data = synthetic_prob(mcmcInfo.seq_length, mcmcInfo.alpha, mcmcInfo.nStates, ...
                          mcmcInfo.nSteps, mcmcInfo.A, mcmcInfo.v, mcmcInfo.noise, mcmcInfo.pi0);                                     
    
    % record full simulation info
    fieldNames = fieldnames(synthetic_data);
    for f = 1:length(fieldNames)
        mcmcInfo.masterSimStruct(i).(fieldNames{f}) = synthetic_data.(fieldNames{f});
    end
end

%% Gibbs MCMC 

% basic inference params 
mcmcInfo.n_mcmc_steps = 100; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.n_chains = 10;

% initialize arrays
mcmcInfo.logL_vec = NaN(mcmcInfo.n_mcmc_steps,1);
mcmcInfo.A_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_mcmc_steps);
mcmcInfo.v_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.n_mcmc_steps);

% specify hyperparameters

% A prior--assume strongly diagonal PDF given short timescale
% take A columns to follow multinomial Dirichlet distribution
A_alpha = ones(mcmcInfo.nStates);%*n_particles*n_traces;
mcmcInfo.A_alpha(eye(mcmcInfo.nStates)==1) = A_alpha(eye(mcmcInfo.nStates)==1)*10; % distribution hyper params

% emission rate priors. For now assume v_i's to be Poisson. Thus, the
% conjugate prior is gamma distributed
v_prior = mcmcInfo.v; % prior on RNAP initiation rates


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo.A_curr = sample_A_dirichlet(A_alpha, zeros(mcmcInfo.nStates));
mcmcInfo.v_curr = mcmcInfo.v;

% calculate pi0 
[V, D] = eig(mcmcInfo.A_curr);
[~, mi] = max(real(diag(D)));
mcmcInfo.pi0_curr = V(:,mi)/sum(V(:,mi));

% initialize chains
mcmcInfo = initialize_chains(mcmcInfo);

% get predicted fluorescence
mcmcInfo = predict_fluo_full(mcmcInfo);


for step = 1:mcmcInfo.n_mcmc_steps
    
    % resample chains
    mcmcInfo = resample_chains(mcmcInfo);
end

%%  