% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
sampling_res = 16.85;
% sampling_res = 2.5;
trueParams = setParamsBasic3state(sampling_res);

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.tres = sampling_res;
mcmcInfo.n_mcmc_steps = 250;
mcmcInfo.n_chains = 25;
mcmcInfo.n_traces = 20;
mcmcInfo.burn_in = 100;
mcmcInfo.seq_length = 100;
mcmcInfo.inferMemory = 0;

mcmcInfo.rateSamplingFlag = 1;
mcmcInfo.adjustSamplingFlag = 0;
mcmcInfo.mhResamplingFlag = 1;

trueParams.discrete_data_flag = 0;
mcmcInfo.resampleTracesFlag = 1;
mcmcInfo.rs_freq = 10;
mcmcInfo.upsample_factor = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = 7;
trueParams.nSteps = nSteps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trueParams.n_traces = mcmcInfo.n_traces;
trueParams.seq_length = mcmcInfo.seq_length;
trueParams = generateSimulatedData(trueParams);

%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%

% add "known" hyperparameters
mcmcInfo.nStates = trueParams.nStates;    
mcmcInfo.alpha_frac = trueParams.alpha_frac;
mcmcInfo.observed_fluo = trueParams.observed_fluo;

if ~mcmcInfo.inferMemory
    mcmcInfo.nSteps = trueParams.nSteps;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
mcmcInfo = setMCMCOptions(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference arrays and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initializeInferenceArrays_v2(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom_v2(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc

%%
close all

figure;
hold on

% plot sampling results
histogram(mcmcInfo.Q_inf_array(1,2,mcmcInfo.burn_in:end,:))
histogram(mcmcInfo.Q_inf_array(2,1,mcmcInfo.burn_in:end,:))
histogram(mcmcInfo.Q_inf_array(3,2,mcmcInfo.burn_in:end,:))
histogram(mcmcInfo.Q_inf_array(2,3,mcmcInfo.burn_in:end,:))


% plot ground truth
yl = get(gca,'ylim');
Q = trueParams.R;
plot([Q(2,1) Q(2,1)], [yl(1) yl(2)],'--k','LineWidth',2);
plot([Q(1,2) Q(1,2)], [yl(1) yl(2)],'--k','LineWidth',2);
plot([Q(2,3) Q(2,3)], [yl(1) yl(2)],'--k','LineWidth',2);
plot([Q(3,2) Q(3,2)], [yl(1) yl(2)],'--k','LineWidth',2);

legend('q12 (true=0.021)','q21 (true=0.02)','q32 (true=0.012)','q23 (true=0.07)')

xlim([0 0.08])