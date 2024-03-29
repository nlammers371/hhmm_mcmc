% Script to validate ns MCMC method for burs parameter inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
if ~exist(DropboxFolder)
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
end    
outPath = [DropboxFolder 'hhmm_MCMC_data\run_mcmc_validation_v3_resampling_test\'];
mkdir(outPath);

iter_size = 50; % parpool deleted and reinitiated every N iterations
%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%
mcmcInfoInit = struct;

% other key hyperparameters
mcmcInfoInit.n_mcmc_steps = 3e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 500;
% mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfoInit.NumWorkers = 12;
mcmcInfoInit.annealingSigmaFlag = 0; % need to implement this
mcmcInfoInit.bootstrapFlag = 0;

% Set the parameter options to explore
seq_length = 100; % length of simulated traces in time steps
inferMemory = 0;
n_chains = 10;
% n_traces_per_chain_vec = [5 10 20 50];
n_traces = 20;
n_resamp_vec = [1 2 5 10];
mem_vec = 7:10;
nStates_vec = [3 2];
simVec = 1:10;
fast3_flag = 0;

sim_struct = struct;
iter = 1;

for m = 1:length(mem_vec)
    for s = 1:length(nStates_vec)

        nStates = nStates_vec(s);
%             n_traces = n_traces_per_chain_vec(n);
        nSteps = mem_vec(m);

        if nStates == 3
            % initialize 3 state system 
            if fast3_flag
                trueParams = setParamsBasic3state_fast;
            else
                trueParams = setParamsBasic3state;
            end
        elseif nStates == 2
            trueParams = setParamsBasic2state;
        else
            error('unsupported number of states')
        end
        trueParams.nSteps = nSteps;
        trueParams.n_traces = n_traces;
        trueParams.seq_length = seq_length;

        % simulate transcription traces
%             sim_struct(iter).trueParams = trueParams;
        sim_struct(iter).trueParams = generateSimulatedData(trueParams);
        sim_struct(iter).nSteps = nSteps;
        sim_struct(iter).n_traces = n_traces;
        sim_struct(iter).nStates = nStates;
        sim_struct(iter).seq_length = seq_length;

        iter = iter + 1;
    end
end  
clear iter

n_step_lookup = [sim_struct.nSteps];
n_state_lookup = [sim_struct.nStates];

% get all possible combinations
elements = {n_resamp_vec mem_vec nStates_vec simVec};
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
        n_reps = combArray(iter,1);
        nSteps = combArray(iter,2);
        nStates = combArray(iter,3);
        rep_num = combArray(iter,4);

        % extract model        
        step_filter = n_step_lookup == nSteps;
        state_filter = n_state_lookup == nStates;

        trueParams = sim_struct(step_filter & state_filter).trueParams;

        % initialize inference results structure
        mcmcInfo = mcmcInfoInit;

        if ~inferMemory
            mcmcInfo.nSteps = nSteps;
        end    

        % add "known" hyperparameters
        mcmcInfo.nStates = trueParams.nStates;  
        mcmcInfo.n_reps = n_reps;  
        mcmcInfo.alpha_frac = trueParams.alpha_frac;
        mcmcInfo.observed_fluo = trueParams.observed_fluo;
        mcmcInfo.n_traces = n_traces;
%         mcmcInfo.n_traces_per_chain = n_traces_per_chain;
        mcmcInfo.seq_length = seq_length;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set MCMC options
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        mcmcInfo = setMCMCOptions(mcmcInfo, n_chains, inferMemory);

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
        saveString = ['K' sprintf('%01d',nStates) '_W' sprintf('%03d',10*round(nSteps,1)) '_nrep' sprintf('%03d',n_reps) '_rep' sprintf('%03d',rep_num)];

        disp('saving...')

        % strip unneccesarry fields    
        mcmcInfo = rmfield(mcmcInfo,{'indArray','trace_logL_array','trace_logL_vec','state_ref','A_curr','v_curr','nStepsCurr','sample_chains_dummy','observed_fluo_dummy','observed_fluo_dummy2','sample_fluo_dummy2'});
        trueParams = rmfield(trueParams,{'masterSimStruct'});
        mcmcInfo.trueParams = trueParams;
        saveFun(mcmcInfo, outPath, saveString)

    end    
end
