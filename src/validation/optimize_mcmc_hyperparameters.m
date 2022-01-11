% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
outPath = [DropboxFolder 'hhmm_MCMC_data\hyperParameterOptimization_2state_v2\'];
mkdir(outPath);

% initialize 3 state system 
mcmcInfoInit = setParamsBasic2state;

% indicate how many replicates of each we want
n_reps = 3;

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfoInit.n_mcmc_steps = 2500; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 500;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step

% characteristics of simulated data
mcmcInfoInit.n_traces = 10;
mcmcInfoInit.seq_length = 120; % length of simulated traces in time steps

%%% Are we doing a consistency test?
mcmcInfoInit.consistencyTestFlag = 0;
mcmcInfoInit.enforceRateConsistency = 0;
mcmcInfoInit = genericInitialization(mcmcInfoInit);

% Set the parameter options to explore
repVec = 1:n_reps;
inferMemoryVec = 0:1;
temperVec = 0:1;
nChainsVec = [10];
nTempsVec = [1 3 5];
nSwaps = [7 14 28];
infoSharing = 0:1;
tempIncrement = [sqrt(2) 2 5];
% get all possible combinations
elements = {inferMemoryVec nChainsVec nTempsVec nSwaps repVec infoSharing tempIncrement temperVec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
combArray = [combCell{:}]; 

% remove irrelevant ones
keepFlags = ~(combArray(:,3)>1 & combArray(:,8)==0) & ~(combArray(:,4)>7 & combArray(:,8)==0) ...
                  & ~(combArray(:,6)==1 & combArray(:,8)==0)  & ~(combArray(:,7)>sqrt(2) & combArray(:,8)==0);
combArray = combArray(keepFlags,:); 

try 
    parpool(24);
catch
    % do nothing
end  

parfor iter = 1:size(combArray,1)
  
    inferMemory = combArray(iter,1)==1;
    n_chains = combArray(iter,2);
    n_temps = combArray(iter,3);
    n_swaps = combArray(iter,4);
    info_sharing = combArray(iter,6);
    temp_increment = combArray(iter,7);
    temperingFlag = combArray(iter,8);
    step_num = combArray(iter,5);
    
    % start timer
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set MCMC options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    mcmcInfo = setMCMCOptions(mcmcInfoInit, n_chains, temperingFlag, n_temps, n_swaps, inferMemory, false,  temp_increment, info_sharing);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo = initializeInferenceArrays(mcmcInfo);
    mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct inference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfo = inferenceWrapper(mcmcInfo);    
    mcmcInfo.duration = toc;
    
    % save results
    saveString = ['nc' sprintf('%03d',n_chains) '_tempering' num2str(temperingFlag) '_ntm' sprintf('%03d',n_temps) '_nsw' sprintf('%03d',n_swaps)...
                '_mem' num2str(inferMemory) '_info' num2str(info_sharing)  '_tempInc' num2str(round(temp_increment,1)) '_rep' sprintf('%03d',step_num)];
    disp('saving...')
    % strip unneccesarry fields    
    mcmcInfo = rmfield(mcmcInfo,{'indArray','trace_logL_array','trace_logL_vec','masterSimStruct','state_ref','A_curr','v_curr','nStepsCurr','sample_chains_dummy','observed_fluo_dummy','observed_fluo_dummy2','sample_fluo_dummy2'});
    saveFun(mcmcInfo, outPath, saveString)
    
end    

