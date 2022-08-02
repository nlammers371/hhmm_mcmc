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
n_mcmc_steps = 250;
n_chains = 25;
n_traces = 20;
seq_length = 100;
inferMemory = 0;
ensembleInferenceFlag = 0;
mcmcInfo.rateSamplingFlag = 1;
trueParams.discrete_data_flag = 0;

% global mcmcInfo
mcmcInfo.n_mcmc_steps = n_mcmc_steps; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 50;
mcmcInfo.resampleTracesFlag = 1;
mcmcInfo.rs_freq = 10;
mcmcInfo.tres = sampling_res;
mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run

% characteristics of simulated data
mcmcInfo.upsample_factor = 2;
mcmcInfo.n_reps = 2;
mcmcInfo.n_traces = n_traces;
mcmcInfo.seq_length = seq_length; % length of simulated traces in time steps

mcmcInfo.mhInferenceFlag = 0;
mcmcInfo.reducedModelFlag = 0;
mcmcInfo.ensembleInferenceFlag = ensembleInferenceFlag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = 7;
trueParams.nSteps = nSteps;

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
mcmcInfo = initializeInferenceArrays_v2(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom_v2(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc