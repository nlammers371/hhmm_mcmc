% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
sampling_res = 16.85;
trueParams = setParamsBasic2state(sampling_res);

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
n_mcmc_steps = 500;
n_chains = 25;
n_traces = 20;
seq_length = 100;
inferMemory = 0;

% global mcmcInfo

mcmcInfo.n_mcmc_steps = n_mcmc_steps; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 150;
mcmcInfo.resampleTracesFlag = 1;
mcmcInfo.rs_freq = 50;
mcmcInfo.tres = sampling_res;
mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run

% characteristics of simulated data
mcmcInfo.upsample_factor = 1;
mcmcInfo.n_reps = 2;
mcmcInfo.n_traces = n_traces;
mcmcInfo.seq_length = seq_length; % length of simulated traces in time steps

mcmcInfo.mhInferenceFlag = 0;
mcmcInfo.reducedModelFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = 7;
trueParams.nSteps = nSteps;
trueParams.discrete_data_flag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trueParams.n_traces = n_traces;
trueParams.seq_length = seq_length;
trueParams = generateSimulatedData(trueParams);

%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%

% add "known" hyperparameters
mcmcInfo.nStates = trueParams.nStates;    
mcmcInfo.alpha_frac = trueParams.alpha_frac;
mcmcInfo.observed_fluo = trueParams.observed_fluo;

% other key hyperparameters
mcmcInfo.n_mcmc_steps = 1e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 500;
mcmcInfo.n_reps = 10; % number of chain state resampling passes per inference step
mcmcInfo.NumWorkers = 4;
mcmcInfo.annealingSigmaFlag = 0; % need to implement this


if ~inferMemory
    mcmcInfo.nSteps = trueParams.nSteps;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
mcmcInfo = setMCMCOptions(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference arrays and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initializeInferenceArrays(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initialize_chains_v3(mcmcInfo);

% get predicted fluorescence
mcmcInfo = predict_fluo_full_v3(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct full inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


mcmcInfo.step = 1;
params_init = [mcmcInfo.k_curr  mcmcInfo.r_curr  mcmcInfo.sigma_curr(1)]';
params_true = [trueParams.kon trueParams.koff log(trueParams.k_corr_factor)...
               0.1 log(1.3) trueParams.sigma]';
loss_init = calculateModelLoss(params_init);
loss_true = calculateModelLoss(params_true);
% mcmcInfo = inferenceWrapper(mcmcInfo);  
loss_fun = @(params) -calculateModelLoss(params);
tic
% options = optimoptions('fmincon','ConstraintTolerance',1);
[params, loss,~,output] = fmincon(loss_fun,params_init,[],[],[],[],[1e-3 1e-3 -0.7 0 -0.7 0],[.1 0.1 0.7 5/mcmcInfo.tres 0.7 5],[]);%,options);
toc