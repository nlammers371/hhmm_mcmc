% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
mcmcInfo = setParamsBasic;

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% characteristics of simulated data
mcmcInfo.n_traces = 12;
mcmcInfo.seq_length = 120; % length of simulated traces in time steps

%%% Are we doing a consistency test?
mcmcInfo.consistencyTestFlag = 0;

%%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
% basic inference params 
mcmcInfo.n_mcmc_steps = 5e3; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 2500;
mcmcInfo.n_chains = 10; % number of parallel MCMC chains to run
mcmcInfo.n_temps_per_chain = 7; % number of rungs in the temperature ladder for each chain
mcmcInfo.n_chains_eff = mcmcInfo.n_temps_per_chain*mcmcInfo.n_chains;
mcmcInfo.n_reps = 1; % number of chain state resampling passes per inference step

% memory parameter
mcmcInfo.trueParams.nSteps = 7; % True parameters
mcmcInfo.nStepsMax = 10; % number of time steps needed to transcribe full gene
mcmcInfo.nStepsCurr = mcmcInfo.trueParams.nSteps; % current guess (can be fractional)
mcmcInfo.alpha = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;

% inference type
mcmcInfo.ensembleInferenceFlag = 0; % perform ensemble inference across parallel chains?

% tempering options
mcmcInfo.temperingFlag = 1; % use parallel tempering?
mcmcInfo.n_rs_per_trace = 7; % number of swap proposals per neighboring trace pair
mcmcInfo.refChainVec = false(1,mcmcInfo.n_chains_eff);
mcmcInfo.refChainVec(1:mcmcInfo.n_temps_per_chain:end) = true; % designate T=1 chains that will actually be used for inference
mcmcInfo.temp_incrememt = 3.5;
exp_vec = repmat(0:mcmcInfo.n_temps_per_chain-1,1,mcmcInfo.n_chains);
mcmcInfo.tempGradVec = mcmcInfo.temp_incrememt.^exp_vec;%logspace(0,log10(mcmcInfo.max_temp),mcmcInfo.n_chains);

% mcmcInfo.tempGradVec(2) = 1;
mcmcInfo.move_flag_array = false(mcmcInfo.n_rs_per_trace,mcmcInfo.n_chains-1,mcmcInfo.n_traces,mcmcInfo.n_mcmc_steps);

% initialize arrays and simulate traces
mcmcInfo = genericInitialization(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref_chain_ids = repelem(1:mcmcInfo.n_chains,mcmcInfo.n_temps_per_chain);
cmap = brewermap(mcmcInfo.n_chains,'Spectral');

close all
n_points = size(mcmcInfo.v_inf_array(:,2,:),1);
x_axis = mcmcInfo.update_increment*(1:n_points);

v3_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,mcmcInfo.v_inf_array(:,3,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,mcmcInfo.v_inf_array(:,3,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.v(3),n_points),'k--','LineWidth',2)
ylabel('state 3 emission rate')
xlabel('MCMC steps')

v2_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,mcmcInfo.v_inf_array(:,2,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,mcmcInfo.v_inf_array(:,2,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.v(2),n_points),'k--','LineWidth',2)
ylabel('state 2 emission rate')
xlabel('MCMC steps')

sigma_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,mcmcInfo.sigma_inf_array(:,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,mcmcInfo.sigma_inf_array(:,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.sigma,n_points),'k--','LineWidth',2)
ylabel('fluorescence noise')
xlabel('MCMC steps')

a21 = NaN(n_points,mcmcInfo.n_chains);
for n = 1:mcmcInfo.n_chains_eff
    a21(:,n) = mcmcInfo.A_inf_array(2,1,:,n);
end
A21_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,a21(:,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,a21(:,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.A(2,1),n_points),'k--','LineWidth',2)

ylabel('pON')
xlabel('MCMC steps')

a12 = NaN(n_points,mcmcInfo.n_chains);
for n = 1:mcmcInfo.n_chains_eff
    a12(:,n) = mcmcInfo.A_inf_array(1,2,:,n);
end

A_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,a12(:,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,a12(:,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.A(1,2),n_points),'k--','LineWidth',2)

ylabel('pOFF')
xlabel('MCMC steps')


% calculate first index to use
first_ind = mcmcInfo.burn_in/mcmcInfo.update_increment + 1;

v3_fig_hist = figure;
v3_data = [];
for i = find(mcmcInfo.refChainVec)
    v3_data = [v3_data mcmcInfo.v_inf_array(first_ind:end,3,i)'];        
end
histogram(v3_data);
% plot(x_axis,repelem(mcmcInfo.trueParams.v(3),n_points),'k--','LineWidth',2)
xlabel('state 3 emission rate')
xlabel('probability')

v2_fig_hist = figure;
v2_data = [];
for i = find(mcmcInfo.refChainVec)
    v2_data = [v2_data mcmcInfo.v_inf_array(first_ind:end,2,i)'];        
end
histogram(v2_data);
% plot(x_axis,repelem(mcmcInfo.trueParams.v(3),n_points),'k--','LineWidth',2)
xlabel('state  emission rate')
xlabel('probability')
%%
sigma_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,mcmcInfo.sigma_inf_array(:,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,mcmcInfo.sigma_inf_array(:,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.sigma,n_points),'k--','LineWidth',2)
ylabel('fluorescence noise')
xlabel('MCMC steps')

a21 = NaN(n_points,mcmcInfo.n_chains);
for n = 1:mcmcInfo.n_chains_eff
    a21(:,n) = mcmcInfo.A_inf_array(2,1,:,n);
end
A21_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,a21(:,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,a21(:,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.A(2,1),n_points),'k--','LineWidth',2)

ylabel('pON')
xlabel('MCMC steps')

a12 = NaN(n_points,mcmcInfo.n_chains);
for n = 1:mcmcInfo.n_chains_eff
    a12(:,n) = mcmcInfo.A_inf_array(1,2,:,n);
end

A_fig = figure;
hold on
for i = 1:mcmcInfo.n_chains_eff
    if mcmcInfo.refChainVec(i)
        plot(x_axis,a12(:,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
    else
        plot(x_axis,a12(:,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
    end
end
plot(x_axis,repelem(mcmcInfo.trueParams.A(1,2),n_points),'k--','LineWidth',2)

ylabel('pOFF')
xlabel('MCMC steps')