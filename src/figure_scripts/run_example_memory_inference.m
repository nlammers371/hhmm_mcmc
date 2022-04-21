% Script to make illustrative inference figures for 2 state case
clear 
close all force

addpath(genpath('../utilities'))

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
rng(958)

n_mcmc_steps = 500;
mcmcInfoInit.n_mcmc_steps = n_mcmc_steps; % number of MCMC steps (need to add convergence criteria)
burn_in = 250;
mcmcInfoInit.burn_in = burn_in;
n_chains = 100;
mcmcInfoInit.n_chains = n_chains; % number of parallel MCMC chains to run
mcmcInfoInit.n_reps = 1; % number of chain state resampling passes per inference step

% characteristics of simulated data
mcmcInfoInit.n_traces = 20;
mcmcInfoInit.seq_length = 120; % length of simulated traces in time steps
mcmcInfoInit.NumWorkers = 5;
% set vector of memories to infer
mem_vec = [4.1 5.7 7.3 8.9 10.4]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initializePool(mcmcInfoInit, 1)

for temp_flag = 1
    mem_array = NaN(n_mcmc_steps+1, n_chains, length(mem_vec));
    logL_array = NaN(length(mem_vec), n_chains);
    parfor n = 1:length(mem_vec)
        % initialize info structure
        trueParams = setParamsBasic3state;
        trueParams.sigma = 1;
        trueParams.discrete_data_flag = false;
        nSteps = mem_vec(n);
        trueParams.nSteps = nSteps;
        mcmcInfo = mcmcInfoInit;
        mcmcInfo.inferNStepsFlag = true;
        mcmcInfo.annealingSigmaFlag = temp_flag==1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulate data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        trueParams.n_traces = mcmcInfo.n_traces;
        trueParams.seq_length = mcmcInfo.seq_length;
        trueParams = generateSimulatedData(trueParams);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % initialize inference arrays and variables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mcmcInfo.nStates = trueParams.nStates;    
        mcmcInfo.alpha_frac = trueParams.alpha_frac;
        mcmcInfo.observed_fluo = trueParams.observed_fluo;

        mcmcInfo = setMCMCOptions(mcmcInfo);

        mcmcInfo = initializeInferenceArrays(mcmcInfo);
        mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

        % mcmcInfo.save_trace_results = false;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % conduct inference
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        mcmcInfo = inferenceWrapper(mcmcInfo);
        toc

        % store results
        mem_array(:,:,n) = mcmcInfo.n_steps_inf_array;
        logL_array(n,:) = mean(mcmcInfo.logL_vec(burn_in+100:end,:),1);
    end
    % save
    mem_data = struct;
    mem_data.logL_array = logL_array;
    mem_data.mem_array = mem_array;
    mem_data.mem_vec = mem_vec;
    disp('saving...')
    save(['mem_data_' num2str(temp_flag) '.mat'],'mem_data')
    disp('done.')
end    