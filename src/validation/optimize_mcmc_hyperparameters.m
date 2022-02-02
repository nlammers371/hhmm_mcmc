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
mcmcInfoInit = setParamsBasic3state;

% indicate how many replicates of each we want
n_sims = 1;

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfoInit.n_mcmc_steps = 1e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 500;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfoInit.swapInc = 10;
mcmcInfoInit.nStepsSwapFlag = 0; % this does not appear to work

% characteristics of simulated data
mcmcInfoInit.n_traces = 2;
mcmcInfoInit.seq_length = 120; % length of simulated traces in time steps

%%% Are we doing a consistency test?
mcmcInfoInit.consistencyTestFlag = 0;
mcmcInfoInit.enforceRateConsistency = 0;
mcmcInfoInit = genericInitialization(mcmcInfoInit);

% Set the parameter options to explore
repVec = 1:n_sims;
inferMemoryVec = 1;
temperVec = 0;
nChainsVec = [10];
nTempsVec = [1];
nSwaps = [14];
infoSharing = 0;
tempIncrement = [1];
% get all possible combinations
elements = {inferMemoryVec nChainsVec nTempsVec nSwaps repVec infoSharing tempIncrement temperVec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
combArray = [combCell{:}]; 

% remove irrelevant ones
% keepFlags = ~(combArray(:,3)>1 & combArray(:,8)==0) & ~(combArray(:,4)>7 & combArray(:,8)==0) ...
%                   & ~(combArray(:,6)==1 & combArray(:,8)==0)  & ~(combArray(:,7)>sqrt(2) & combArray(:,8)==0);
% combArray = combArray(keepFlags,:); 

try 
    parpool(4);
catch
    % do nothing
end  

for iter = 1:size(combArray,1)
  
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
    mcmcInfoTemp = setMCMCOptions(mcmcInfoInit, n_chains, temperingFlag, n_temps, n_swaps, inferMemory, false,  temp_increment, info_sharing);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo = mcmcInfoTemp;
%     mcmcInfoTemp.n_mcmc_steps = mcmcInfoInit.burn_in;
%     mcmcInfoTemp = initializeInferenceArrays(mcmcInfoTemp);
%     mcmcInfoTemp = initializeVariablesBasicRandom(mcmcInfoTemp);
%     tic
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % conduct initial inference
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%     
%     mcmcInfoTemp = inferenceWrapper(mcmcInfoTemp);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct full inference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfo = initializeInferenceArrays(mcmcInfo);
%     ids_to_use = randsample(1:mcmcInfoTemp.n_chains,mcmcInfoTemp.n_chains,true,exp(sum(mcmcInfoTemp.trace_logL_vec,3)));    
%     mcmcInfo = initializeVariablesWithPriors(mcmcInfo,mcmcInfoTemp,ids_to_use);
    mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);
    mcmcInfo = inferenceWrapper(mcmcInfo);  
    
    mcmcInfo.duration = toc;
    
    % save results
%     saveString = ['nc' sprintf('%03d',n_chains) '_tempering' num2str(temperingFlag) '_ntm' sprintf('%03d',n_temps) '_nsw' sprintf('%03d',n_swaps)...
%                 '_mem' num2str(inferMemory) '_info' num2str(info_sharing)  '_tempInc' num2str(round(temp_increment,1)) '_rep' sprintf('%03d',step_num)];
    saveString = ['nc' sprintf('%03d',n_chains) '_tempering' num2str(temperingFlag) '_ntm' sprintf('%03d',n_temps) '_nsw' sprintf('%03d',n_swaps)...
                '_mem' num2str(inferMemory) '_tempInc' num2str(round(10*temp_increment,0)) '_rep' sprintf('%03d',step_num)];

    disp('saving...')
    % strip unneccesarry fields    
    mcmcInfo = rmfield(mcmcInfo,{'indArray','trace_logL_array','trace_logL_vec','masterSimStruct','state_ref','A_curr','v_curr','nStepsCurr','sample_chains_dummy','observed_fluo_dummy','observed_fluo_dummy2','sample_fluo_dummy2'});
    saveFun(mcmcInfo, outPath, saveString)
    
end    

