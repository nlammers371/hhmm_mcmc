% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = setParamsBasic;

%%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_traces = 11;
mcmcInfo.n_mcmc_steps = 200; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.update_increment = 10; % sets how often parameter values are recorded in inference arrays
mcmcInfo.n_chains = 20;

% initialize arrays and simulate traces
mcmcInfo = genericInitialization(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = intitializeVariablesBasicRandom(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc
