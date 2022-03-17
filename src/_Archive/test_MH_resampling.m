% Script to test that model returns correct parameters when initiated with
% correct parameters
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = setParamsBasic;

% set key system parameters
mcmcInfo.n_traces = 5;
mcmcInfo.nSteps = 7;
mcmcInfo.seq_length = 120*60/mcmcInfo.tres;
mcmcInfo.testResampling = 1;
mcmcInfo.n_reps = 1;
mcmcInfo.MHResampling = 1;

%%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_mcmc_steps = 50; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.update_increment = 1; % sets how often parameter values are recorded in inference arrays
mcmcInfo.n_chains = 10;

% initialize arrays and simulate traces
mcmcInfo = genericInitialization(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mcmcInfo = intitializeVariablesTrue(mcmcInfo);
mcmcInfo = intitializeVariablesBasicRandom(mcmcInfo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc
