% Script to test number of resampling steps needed to obtain optimal fit
% given correct parameters
clear 
close all force

addpath(genpath('utilities'))
%%
% initialize info structure
mcmcInfo = setParamsBasic;

% set key system parameters
mcmcInfo.n_traces = 3;
mcmcInfo.nSteps = 6;
mcmcInfo.seq_length = 120*60/mcmcInfo.tres;
mcmcInfo.testResampling = 1;
mcmcInfo.n_reps = 25;

%%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_mcmc_steps = 2; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.update_increment = 10; % sets how often parameter values are recorded in inference arrays
mcmcInfo.n_chains = 5;

% initialize arrays and simulate traces
mcmcInfo = genericInitialization(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = intitializeVariablesTrue(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc
