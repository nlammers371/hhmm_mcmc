% Script to validate ns MCMC method for burs parameter inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
if ~exist(DropboxFolder)
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
end    
outPath = [DropboxFolder 'hhmm_MCMC_data\mcmc_testomg\'];
mkdir(outPath);

iter_size = 50; % parpool deleted and reinitiated every N iterations
%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%

mcmcInfoInit = struct;

% other key hyperparameters
mcmcInfoInit.mhResamplingFlag = 1;
mcmcInfoInit.n_mcmc_steps = 500; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 2e3;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfoInit.NumWorkers = 25;
mcmcInfoInit.ensembleInferenceFlag = 1;

% Set the parameter options to explore
seq_length = 100; % length of simulated traces in time steps
mcmcInfoInit.n_chains = 25;
n_traces = 20;
nSteps = 7;
nStates = 3;
% n_reps = 5;
master_struct = struct;

% initialize parallel pool
initializePool(mcmcInfoInit, 0)

for iter = 1
  
    % generate model
    if nStates == 3
        % initialize 3 state system   
        trueParams = setParamsBasic3state;        
    elseif nStates == 2
        trueParams = setParamsBasic2state;
    else
        error('unsupported number of states')
    end
    trueParams.nSteps = nSteps;
    trueParams.n_traces = n_traces;
    trueParams.seq_length = seq_length;

    % simulate traces
    trueParams = generateSimulatedData(trueParams);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set MCMC options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % initialize inference results structure
    mcmcInfo = mcmcInfoInit;

    % add "known" hyperparameters                
    mcmcInfo.nStates = trueParams.nStates;    
    mcmcInfo.alpha_frac = trueParams.alpha_frac;
    mcmcInfo.observed_fluo = trueParams.observed_fluo;
    mcmcInfo.n_traces = n_traces;
    mcmcInfo.n_traces_per_chain = n_traces;
    mcmcInfo.seq_length = seq_length;

    % check and implement options
    mcmcInfo = setMCMCOptions(mcmcInfo);       

    if ~mcmcInfo.inferNStepsFlag
        mcmcInfo.nSteps = nSteps;
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo = initializeInferenceArrays(mcmcInfo);
    mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct full inference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    tic    
    mcmcInfo = inferenceWrapper(mcmcInfo);      
    mcmcInfo.duration = toc;   
    mcmcInfo.trueParams = trueParams;
    master_struct(iter).mcmcInfo = mcmcInfo;
    master_struct(iter).logL = mean(mcmcInfo.logL_vec(end,:),1);

end    

