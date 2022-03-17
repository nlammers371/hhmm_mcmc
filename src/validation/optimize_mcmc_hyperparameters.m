% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
if ~exist(DropboxFolder)
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
end    
outPath = [DropboxFolder 'hhmm_MCMC_data\hyperParameterOptimization_3state\'];
mkdir(outPath);

% initialize 3 state system 
trueParams = setParamsBasic3state;

% characteristics of simulated data
trueParams.n_traces = 10;
trueParams.seq_length = 120; % length of simulated traces in time steps

% simulate transcription traces
trueParams = generateSimulatedData(trueParams);

% indicate how many replicates of each we want
n_sims = 1;

%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%
mcmcInfoInit = struct;

% add "known" hyperparameters
mcmcInfoInit.nStates = trueParams.nStates;
mcmcInfoInit.alpha_frac = trueParams.alpha_frac;
mcmcInfoInit.observed_fluo = trueParams.observed_fluo;
mcmcInfoInit.n_traces = size(mcmcInfoInit.observed_fluo,2);
mcmcInfoInit.seq_length = size(mcmcInfoInit.observed_fluo,1);

% other key hyperparameters
mcmcInfoInit.n_mcmc_steps = 1e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 500;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfoInit.nStepsSwapFlag = 0; % this does not appear to work
mcmcInfoInit.NumWorkers = 4;
mcmcInfoInit.annealingSigmaFlag = 0; % need to implement this

% Set the parameter options to explore
repVec = 1:n_sims;
inferMemoryVec = 0;
nChainsVec = [10];


if ~inferMemoryVec
    mcmcInfoInit.nSteps = trueParams.nSteps;
end    
% get all possible combinations
elements = {inferMemoryVec nChainsVec repVec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
combArray = [combCell{:}]; 

% initialize parallel pool
initializePool(mcmcInfoInit)

for iter = 1:size(combArray,1)
  
    inferMemory = combArray(iter,1)==1;
    n_chains = combArray(iter,2);
%     n_temps = combArray(iter,3);
%     n_swaps = combArray(iter,4);
%     info_sharing = combArray(iter,6);
%     temp_increment = combArray(iter,7);
%     temperingFlag = combArray(iter,8);
    step_num = combArray(iter,end);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set MCMC options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    mcmcInfoTemp = setMCMCOptions(mcmcInfoInit, n_chains, inferMemory);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo = mcmcInfoTemp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct full inference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    tic
    mcmcInfo = initializeInferenceArrays(mcmcInfo);
    mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);
    mcmcInfo = inferenceWrapper(mcmcInfo);  
    
    mcmcInfo.duration = toc;
    
    % save results
    saveString = ['nc' sprintf('%03d',n_chains) '_tempering' num2str(temperingFlag) '_ntm' sprintf('%03d',n_temps) '_nsw' sprintf('%03d',n_swaps)...
                '_mem' num2str(inferMemory) '_tempInc' num2str(round(10*temp_increment,0)) '_rep' sprintf('%03d',step_num)];

    disp('saving...')
    % strip unneccesarry fields    
    mcmcInfo = rmfield(mcmcInfo,{'indArray','trace_logL_array','trace_logL_vec','masterSimStruct','state_ref','A_curr','v_curr','nStepsCurr','sample_chains_dummy','observed_fluo_dummy','observed_fluo_dummy2','sample_fluo_dummy2'});
    saveFun(mcmcInfo, outPath, saveString)
    
end    

