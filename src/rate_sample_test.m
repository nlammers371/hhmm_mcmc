% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
sampling_res = 16.85;
% sampling_res = 2.5;
trueParams = setParamsBasic2state(sampling_res);

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.tres = sampling_res;
mcmcInfo.n_mcmc_steps = 150;
mcmcInfo.n_chains = 25;
mcmcInfo.n_traces = 20;
mcmcInfo.burn_in = 50;
mcmcInfo.seq_length = 100;
mcmcInfo.inferMemory = 0;
mcmcInfo.rateSamplingFlag = 0;
trueParams.discrete_data_flag = 0;
mcmcInfo.resampleTracesFlag = 1;
mcmcInfo.rs_freq = 10;
mcmcInfo.upsample_factor = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = 7;
trueParams.nSteps = nSteps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trueParams.n_traces = mcmcInfo.n_traces;
trueParams.seq_length = mcmcInfo.seq_length;
trueParams = generateSimulatedData(trueParams);

%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%

% add "known" hyperparameters
mcmcInfo.nStates = trueParams.nStates;    
mcmcInfo.alpha_frac = trueParams.alpha_frac;
mcmcInfo.observed_fluo = trueParams.observed_fluo;

if ~mcmcInfo.inferMemory
    mcmcInfo.nSteps = trueParams.nSteps;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
mcmcInfo = setMCMCOptions(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference arrays and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initializeInferenceArrays_v2(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom_v2(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc