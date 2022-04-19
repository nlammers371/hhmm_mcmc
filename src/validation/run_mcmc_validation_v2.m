% Script to validate ns MCMC method for burs parameter inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
if ~exist(DropboxFolder)
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
end    
outPath = [DropboxFolder 'hhmm_MCMC_data\mcmc_validation_basic_v2\'];
mkdir(outPath);

iter_size = 50; % parpool deleted and reinitiated every N iterations
%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%
mcmcInfoInit = struct;

% other key hyperparameters
mcmcInfoInit.mhResamplingFlag = 0;
mcmcInfoInit.n_mcmc_steps = 5e0; % number of MCMC steps (need to add convergence criteria)
mcmcInfoInit.burn_in = 2e3;
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfoInit.NumWorkers = 25;

% Set the parameter options to explore
seq_length = 120; % length of simulated traces in time steps
total_chains_per_sim = 125;
n_chains = 5;
n_mcmc_runs = total_chains_per_sim/n_chains;
n_traces = 25;
nSteps = 7;
discrete_data_vec = [0 1];
n_rep_vec = [1 2];
nStates_vec = [3 2];
simVec = 1:10;

% get all possible combinations
elements = {n_rep_vec discrete_data_vec nStates_vec simVec};
combCell = cell(1, numel(elements));
[combCell{:}] = ndgrid(elements{:});
combCell = cellfun(@(x) x(:), combCell,'uniformoutput',false); %there may be a better way to do this
combArray = [combCell{:}]; 

n_blocks = ceil(size(combArray,1)/iter_size);

wb = waitbar(0,'Running MCMC validations');
 
for iter = 1%:size(combArray,1)
    waitbar(iter/size(combArray,1),wb);
    % extract sim characteristics 
    n_reps = combArray(iter,1);
    discrete_data_flag = combArray(iter,2);
    nStates = combArray(iter,3);
    sim_num = combArray(iter,4);
    
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
    trueParams.discrete_data_flag = discrete_data_flag;
    
    % simulate traces
    trueParams = generateSimulatedData(trueParams);
        
    % initialize parallel pool
    initializePool(mcmcInfoInit, 1)
    mcmcInfoFull = struct;
    
   
    parfor n = 1:n_mcmc_runs
                                                      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set MCMC options
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
        % initialize inference results structure
        mcmcInfoTemp = mcmcInfoInit;
        
        % add "known" hyperparameters                
        mcmcInfoTemp.nStates = trueParams.nStates;    
        mcmcInfoTemp.alpha_frac = trueParams.alpha_frac;
        mcmcInfoTemp.observed_fluo = trueParams.observed_fluo;
        mcmcInfoTemp.n_traces = n_traces;
        mcmcInfoTemp.n_traces_per_chain = n_traces;
        mcmcInfoTemp.seq_length = seq_length;
        mcmcInfoTemp.n_reps = n_reps;
        
        % check and implement options
        mcmcInfoTemp = setMCMCOptions(mcmcInfoTemp);       
        
        if ~mcmcInfoTemp.inferNStepsFlag
            mcmcInfoTemp.nSteps = nSteps;
        end    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initialize inference arrays and variables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mcmcInfoTemp = initializeInferenceArrays(mcmcInfoTemp);
        mcmcInfoTemp = initializeVariablesBasicRandom(mcmcInfoTemp);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conduct full inference
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        tic    
        mcmcInfoTemp = inferenceWrapper(mcmcInfoTemp);      
        mcmcInfoTemp.duration = toc;

        fields = fieldnames(mcmcInfoTemp);
        for k = 1:numel(fields)
            aField = fields{k};
            mcmcInfoFull(n).(aField)= mcmcInfoTemp.(aField);        
        end
    end    
    
    % save results
    saveString = ['K' sprintf('%01d',nStates) '_W' sprintf('%03d',10*round(nSteps,1)) '_nt' sprintf('%03d',n_traces) ...
                  '_nrep' sprintf('%1d',n_reps) '_disc' sprintf('%1d',discrete_data_flag)  '_sim' sprintf('%03d',sim_num)];

    disp('saving...')
    fields = fieldnames(mcmcInfoFull);
    mcmcInfo = struct;
    for k = 1:numel(fields)
      aField = fields{k}; % EDIT: changed to {}
      dim = length(size(mcmcInfoFull(1).(aField)));
%       if dim == 3
%           cat_dim = 3;
%       elseif dim == 4
%           cat_dim = 4;
%       else
%           cat_dim = 2;
%       end
      mcmcInfo.(aField) = cat(dim, mcmcInfoFull.(aField));
    end
    mcmcInfo = [mcmcInfoFull.mcmcInfo];
    % strip unneccesarry fields            
    mcmcInfo = rmfield(mcmcInfo,{'trace_logL_array','trace_logL_vec','state_ref','A_curr','v_curr','nStepsCurr'});
    trueParams = rmfield(trueParams,{'masterSimStruct'});
    mcmcInfo.trueParams = trueParams;
    saveFun(mcmcInfo, outPath, saveString)
end
delete(wb)