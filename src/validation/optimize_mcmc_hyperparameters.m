% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
outPath = '../../out/hyperParameterOptimization/';
mkdir(outPath);

% initialize 3 state system 
mcmcInfoInit = setParamsBasic3state;

% indicate how many replicates of each we want
n_reps = 10;

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfoInit.n_mcmc_steps = 5e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 500;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step

% characteristics of simulated data
mcmcInfoInit.n_traces = 20;
mcmcInfoInit.seq_length = 120; % length of simulated traces in time steps

%%% Are we doing a consistency test?
mcmcInfoInit.consistencyTestFlag = 0;
mcmcInfoInit.enforceRateConsistency = 0;
mcmcInfoInit = genericInitialization(mcmcInfoInit);

% Set the parameter options to explore
repVec = 1:n_reps;
inferMemoryVec = 0:1;
nChainsVec = [1 5:10:45];
nTempsVec = [1 3 5 7];
nSwaps = [1 7 14 28];
% get all possible combinations
elements = {inferMemoryVec nChainsVec nTempsVec nSwaps repVec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
combArray = [combCell{:}]; 

% remove irrelevant ones
keepFlags = combArray(:,3)>1 | combArray(:,4)==1;
combArray = combArray(keepFlags,:); 

% try 
%     parpool(24);
% catch
%     % do nothing
% end  

for iter = 1:size(combArray,1)
  
    inferMemory = combArray(iter,1)==1;
    n_chains = combArray(iter,2);
    n_temps = combArray(iter,3);
    temperingFlag = n_temps>1;
    n_swaps = combArray(iter,end-1);
    step_num = combArray(iter,end);
    
    % start timer
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set MCMC options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfo = setMCMCOptions(mcmcInfoInit, n_chains, temperingFlag, n_temps, n_swaps, inferMemory);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo = initializeInferenceArrays(mcmcInfo);
    mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct inference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfo = inferenceWrapper(mcmcInfo);
    mcmcInfo.duration = toc;
    
    % save results
    saveString = ['nc' num2str(n_chains) '_ntm' num2str(n_temps) '_nsw' num2str(n_swaps) '_mem' num2str(inferMemory) '_rep' num2str(step_num)];
    saveFun(mcmcInfo, outPath, saveString)
        
end    
