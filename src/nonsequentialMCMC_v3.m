% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = setParamsBasic;

%%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_traces = 5;
mcmcInfo.n_mcmc_steps = 200; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.update_increment = 10; % sets how often parameter values are recorded in inference arrays
mcmcInfo.n_chains = 20;
% mcmcInfo.nSteps = 10;
mcmcInfo.seq_length = 120;
mcmcInfo.n_reps = 1;

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