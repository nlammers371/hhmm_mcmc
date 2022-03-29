% Script to validate ns MCMC method for burs parameter inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
if ~exist(DropboxFolder)
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
end    
outPath = [DropboxFolder 'hhmm_MCMC_data\mcmc_validation_basic\'];
mkdir(outPath);

iter_size = 50; % parpool deleted and reinitiated every N iterations
%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%
mcmcInfoInit = struct;

% other key hyperparameters
mcmcInfoInit.n_mcmc_steps = 5e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 2e3;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfoInit.NumWorkers = 25;

% Set the parameter options to explore
seq_length = 100; % length of simulated traces in time steps
inferMemory = 0;
mcmcInfo.n_chains = 25;
n_traces_vec = [10 20 50];
mem_vec = 6:10;
nStates_vec = [3 2];
simVec = 1:10;

% get all possible combinations
elements = {n_traces_vec mem_vec nStates_vec simVec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
combArray = [combCell{:}]; 

n_blocks = ceil(size(combArray,1)/iter_size);

for n = 1:n_blocks
    iter_min = (n-1)*iter_size+1;
    iter_max = min(n*iter_size,size(combArray,1));
    
    % initialize parallel pool
    initializePool(mcmcInfoInit, 1)
    
    parfor iter = iter_min:iter_max

        % extract sim characteristics
        n_traces = combArray(iter,1);
        nSteps = combArray(iter,2);
        nStates = combArray(iter,3);
        rep_num = combArray(iter,4);

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

        % save results
        saveString = ['K' sprintf('%01d',nStates) '_W' sprintf('%03d',10*round(nSteps,1)) '_nt' sprintf('%03d',n_traces) '_rep' sprintf('%03d',rep_num)];

        disp('saving...')

        % strip unneccesarry fields            
        mcmcInfo = rmfield(mcmcInfo,{'trace_logL_array','trace_logL_vec','state_ref','A_curr','v_curr','nStepsCurr'});
        trueParams = rmfield(trueParams,{'masterSimStruct'});
        mcmcInfo.trueParams = trueParams;
        saveFun(mcmcInfo, outPath, saveString)

    end    
end
