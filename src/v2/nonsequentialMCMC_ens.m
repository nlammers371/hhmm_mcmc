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
mcmcInfo.v = [.05, 2, 4]';
mcmcInfo.seq_length = 120*60/mcmcInfo.tres;

mcmcInfo.nSteps = 3;
[V, D] = eig(mcmcInfo.A);
[~, mi] = max(real(diag(D)));
mcmcInfo.pi0 = V(:,mi)/sum(V(:,mi));

mcmcInfo.sigma = 0.5;
mcmcInfo.alpha = 1;
mcmcInfo.n_traces = 10;
mcmcInfo.eps = 1e-2; % NL: note that this is not currently used

%%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_mcmc_steps = 250; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.update_increment = 10; % sets how often parameter values are recorded in inference arrays
mcmcInfo.n_chains = 10;

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
n_updates = mcmcInfo.n_mcmc_steps/mcmcInfo.update_increment + 1;
mcmcInfo.logL_vec = NaN(n_updates,1);
mcmcInfo.A_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.nStates,n_updates);
mcmcInfo.v_inf_array = NaN(n_updates,mcmcInfo.nStates);
mcmcInfo.sigma_inf_array = NaN(n_updates,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A prior--assume strongly diagonal PDF given short timescale
% take A columns to follow multinomial Dirichlet distribution
mcmcInfo.A_alpha = ones(mcmcInfo.nStates);%*n_particles*n_traces;
mcmcInfo.A_alpha(eye(mcmcInfo.nStates)==1) = mcmcInfo.A_alpha(eye(mcmcInfo.nStates)==1) + rand(mcmcInfo.nStates,1)*10; % distribution hyper params
mcmcInfo.A_curr = sample_A_dirichlet_par(mcmcInfo.A_alpha, zeros(mcmcInfo.nStates), 1);
mcmcInfo.A_inf_array(:,:,1) = mcmcInfo.A_curr;

% calculate pi0 
mcmcInfo.pi0_curr = NaN(mcmcInfo.n_chains,mcmcInfo.nStates);
[V, D] = eig(mcmcInfo.A_curr);
[~, mi] = max(real(diag(D)));
mcmcInfo.pi0_curr = V(:,mi)/sum(V(:,mi));
mcmcInfo.pi0_inf_array(1,:) = mcmcInfo.pi0_curr;

% initialize sigma as inverse gamma (see: http://ljwolf.org/teaching/gibbs.html)
fluo_vec = mcmcInfo.observed_fluo(:);
f_factor = 0.1*mean(fluo_vec);
mcmcInfo.sigma_curr = trandn(-1,Inf)*f_factor/2 + f_factor;%sqrt(1./gamrnd(100*mcmcInfo.seq_length*mcmcInfo.n_traces/2,1./(fluo_vec'*fluo_vec)));
mcmcInfo.sigma_inf_array(1) = mcmcInfo.sigma_curr;

% initialize v
v2 = prctile(fluo_vec,99) / mcmcInfo.nSteps;%mean(fluo_vec)/sum(mcmcInfo.coeff_MS2)/(mcmcInfo.pi0_curr(2)+2*mcmcInfo.pi0_curr(3));
mcmcInfo.v_curr = [0 v2 2*v2]' + rand(mcmcInfo.nStates,1) - .5;
mcmcInfo.v_inf_array(1,:) = mcmcInfo.v_curr;

% initialize chains
mcmcInfo = initialize_chains_ens(mcmcInfo);

% get predicted fluorescence
mcmcInfo = predict_fluo_full(mcmcInfo);

wb = waitbar(0,'conducting MCMC inference...');
tic
for step = 2:mcmcInfo.n_mcmc_steps %mcmcInfo.n_mcmc_steps    
    waitbar(step/mcmcInfo.n_mcmc_steps,wb);
    
    mcmcInfo.step = step;
    
    % resample chains    
    mcmcInfo = resample_chains_ens(mcmcInfo);              
    
    % get empirical transition and occupancy counts    
    mcmcInfo = get_empirical_counts_ens(mcmcInfo);
   
    % use Gibbs sampling to update hyperparameters      
    mcmcInfo = update_hmm_parameters_ens(mcmcInfo);    
 
end
toc
disp('done')
delete(wb);