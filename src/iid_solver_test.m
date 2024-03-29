% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('../utilities'))

% initialize 3 state system 
trueParams = setParamsBasic3state;

% characteristics of simulated data
trueParams.n_traces = 20;
trueParams.seq_length = 120; % length of simulated traces in time steps

% simulate transcription traces
trueParams = generateSimulatedData(trueParams);

% indicate how many replicates of each we want
n_sims = 1;

%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%
mcmcInfo = struct;

% add "known" hyperparameters
mcmcInfo.nStates = trueParams.nStates;
mcmcInfo.alpha_frac = trueParams.alpha_frac;
mcmcInfo.observed_fluo = trueParams.observed_fluo;
mcmcInfo.n_traces = size(mcmcInfo.observed_fluo,2);
mcmcInfo.seq_length = size(mcmcInfo.observed_fluo,1);

% other key hyperparameters
mcmcInfo.n_mcmc_steps = 3e2; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 125;
mcmcInfo.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfo.nStepsSwapFlag = 0; % this does not appear to work
mcmcInfo.NumWorkers = 4;
mcmcInfo.annealingSigmaFlag = 0; % need to implement this

% Set the parameter options to explore
repVec = 1:n_sims;
inferMemory = 0;
n_chains = 25;

if ~inferMemory
    mcmcInfo.nSteps = trueParams.nSteps;
end    

% initialize parallel pool
% initializePool(mcmcInfo)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
mcmcInfo = setMCMCOptions(mcmcInfo, n_chains, inferMemory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference arrays and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initializeInferenceArrays(mcmcInfo);
% mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);
mcmcInfo = initializeVariables_iid(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = inferenceWrapper_iid(mcmcInfo);