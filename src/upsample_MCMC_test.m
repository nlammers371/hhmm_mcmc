% Script to validate ns MCMC method for burs parameter inference
clear 
close all force

addpath(genpath('./utilities'))

%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%
mcmcInfo = struct;

% other key hyperparameters
mcmcInfo.n_mcmc_steps = 5e2; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 100;
mcmcInfo.n_reps = 1; % number of chain state resampling passes per inference step
% mcmcInfoInit.NumWorkers = 24;
mcmcInfo.annealingSigmaFlag = 0; % need to implement this

% initialize model
trueParams = setParamsBasic3state;

% set system characteristics
seq_length_obs = 100; % length of simulated traces in time steps
inferMemory = 0;
n_chains = 10;
n_traces = 20;
trueParams.nSteps_obs = 7.1;
trueParams.upsample_factor = 1;
trueParams.nSteps = trueParams.nSteps*trueParams.upsample_factor;
mcmcInfo.tres_obs = trueParams.upsample_factor*trueParams.tres; % "true" resolution of system

trueParams.n_traces = n_traces;
trueParams.seq_length = trueParams.upsample_factor*seq_length_obs;
            
% simulate transcription traces
trueParams = generateSimulatedData(trueParams);

if ~inferMemory
    mcmcInfo.nSteps = trueParams.nSteps;
    mcmcInfo.nSteps_obs = mcmcInfo.nSteps/trueParams.upsample_factor;
end    

% add "known" hyperparameters
mcmcInfo.nStates = trueParams.nStates;    
mcmcInfo.alpha_frac = trueParams.alpha_frac;
mcmcInfo.observed_fluo = trueParams.observed_fluo;
mcmcInfo.n_traces = n_traces;
mcmcInfo.seq_length_obs = seq_length_obs;
mcmcInfo.seq_length = trueParams.seq_length;
mcmcInfo.upsample_factor = trueParams.upsample_factor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
mcmcInfo = setMCMCOptions(mcmcInfo, n_chains, inferMemory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference arrays and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initializeInferenceArrays(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
tic    
mcmcInfo = inferenceWrapper(mcmcInfo);      
mcmcInfo.duration = toc;

