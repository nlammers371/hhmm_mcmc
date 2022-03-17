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
mcmcInfoInit.n_mcmc_steps = 2.5e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 250;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfoInit.swapInc = 10; % NL: not used atm

bootstrapFlag = 0;
annealingFlag = 0;
annealingSigmaFlag = 1;

% characteristics of simulated data
mcmcInfoInit.n_traces = 10;
mcmcInfoInit.seq_length = 120; % length of simulated traces in time steps

%%% Are we doing a consistency test?
mcmcInfoInit.consistencyTestFlag = 0; % NL: not used currently
mcmcInfoInit = simulateRawData(mcmcInfoInit);

% Set the parameter options to explore
repVec = 1:n_sims;
inferMemoryVec = 1;
temperVec = 0;
nChainsVec = [25];
nTempsVec = [1]; 
nSwaps = [14]; % NL: not used currently
tempIncrement = [1]; 

% get all possible combinations
elements = {inferMemoryVec nChainsVec nTempsVec nSwaps repVec 0 tempIncrement temperVec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
combArray = [combCell{:}]; 

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set MCMC options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    mcmcInfoAnneal = setMCMCOptions(mcmcInfoInit, n_chains, temperingFlag, n_temps, n_swaps, inferMemory, false,...
                              temp_increment, info_sharing, bootstrapFlag, annealingFlag, annealingSigmaFlag);                                                       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    mcmcInfoAnneal = initializeInferenceArrays(mcmcInfoAnneal);
    mcmcInfoAnneal = initializeVariablesBasicRandom(mcmcInfoAnneal);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct full inference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfoAnneal = inferenceWrapper_v2(mcmcInfoAnneal);  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    mcmcInfo = setMCMCOptions(mcmcInfoInit, n_chains, temperingFlag, n_temps, n_swaps, inferMemory, false,...
                              temp_increment, info_sharing, bootstrapFlag, 0, 0);                                                       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % generate prior
    trunc_array = mcmcInfoAnneal.n_steps_inf_array(mcmcInfoAnneal.burn_in:end,:);
    mu = mean(trunc_array(:));
    sigma = std(trunc_array(:));
    prior_dist = makedist('Normal',mu,sigma);
    mcmcInfo.nSteps_prior = truncate(prior_dist,mcmcInfo.nStepsMin,mcmcInfo.nStepsMax);
    
    mcmcInfo = initializeInferenceArrays(mcmcInfo);
    mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);
    
    mcmcInfo = inferenceWrapper_v2(mcmcInfo);      
    
    % save results
    saveString = ['nc' sprintf('%03d',n_chains) '_tempering' num2str(temperingFlag) '_ntm' sprintf('%03d',n_temps) '_nsw' sprintf('%03d',n_swaps)...
                '_mem' num2str(inferMemory) '_tempInc' num2str(round(10*temp_increment,0)) '_rep' sprintf('%03d',step_num)];

    disp('saving...')
    
    % strip unneccesarry fields    
    mcmcInfoAnneal = rmfield(mcmcInfoAnneal,{'indArray','trace_logL_array','trace_logL_vec','masterSimStruct','state_ref','A_curr','v_curr','nStepsCurr','sample_chains_dummy','observed_fluo_dummy','observed_fluo_dummy2','sample_fluo_dummy2'});
    saveFun(mcmcInfoAnneal, outPath, saveString)
    
end    

