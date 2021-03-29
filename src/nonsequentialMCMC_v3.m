% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = setParamsBasic;

%%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% characteristics of simulated data
mcmcInfo.n_traces = 10;
mcmcInfo.seq_length = 120; % length of simulated traces in time steps

% basic inference params 
mcmcInfo.n_mcmc_steps = 3e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.n_chains = 20; % number of parallel MCMC chains to run
mcmcInfo.n_reps = 1; % number of chain state resampling passes per inference step
mcmcInfo.nSteps = 10; % Number of time steps needed to transcribe full gene

% inference type
mcmcInfo.ensembleInferenceFlag = 0; % perform ensemble inference across parallel chains?
mcmcInfo.temperingFlag = 1; % use parallel tempering?
mcmcInfo.n_rs_per_trace = 5;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
n_points = size(mcmcInfo.v_inf_array(:,2,:),1);
x_axis = 10*(1:n_points);
v3_fig = figure;
hold on
plot(x_axis,permute(mcmcInfo.v_inf_array(:,3,:),[1 3 2]),'-o')
plot(x_axis,repelem(mcmcInfo.trueParams.v(3),n_points),'k--','LineWidth',2)
ylabel('state 3 emission rate')
xlabel('MCMC steps')

v2_fig = figure;
hold on
plot(x_axis,permute(mcmcInfo.v_inf_array(:,2,:),[1 3 2]),'-o')
plot(x_axis,repelem(mcmcInfo.trueParams.v(2),n_points),'k--','LineWidth',2)
ylabel('state 2 emission rate')
xlabel('MCMC steps')

sigma_fig = figure;
hold on
plot(x_axis,mcmcInfo.sigma_inf_array,'-o')
plot(x_axis,repelem(mcmcInfo.trueParams.sigma,n_points),'k--','LineWidth',2)
ylabel('fluorescence noise')
xlabel('MCMC steps')

a21 = NaN(n_points,mcmcInfo.n_chains);
for n = 1:mcmcInfo.n_chains
    a21(:,n) = mcmcInfo.A_inf_array(2,1,:,n);
end
A21_fig = figure;
hold on
plot(x_axis,a21,'-o')
plot(x_axis,repelem(mcmcInfo.trueParams.A(2,1),n_points),'k--','LineWidth',2)

ylabel('pON')
xlabel('MCMC steps')

a12 = NaN(n_points,mcmcInfo.n_chains);
for n = 1:mcmcInfo.n_chains
    a12(:,n) = mcmcInfo.A_inf_array(1,2,:,n);
end

A_fig = figure;
hold on
plot(x_axis,a12,'-o')
plot(x_axis,repelem(mcmcInfo.trueParams.A(1,2),n_points),'k--','LineWidth',2)

ylabel('pOFF')
xlabel('MCMC steps')