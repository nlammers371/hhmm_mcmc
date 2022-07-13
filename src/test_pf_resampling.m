% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
trueParams = setParamsBasic2state;

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
n_mcmc_steps = 100;
mcmcInfo.n_mcmc_steps = n_mcmc_steps; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 100;
n_chains = 101;
mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run

% characteristics of simulated data
mcmcInfo.n_reps = 2;
mcmcInfo.n_traces = 20;
mcmcInfo.n_traces_per_chain = 10;
mcmcInfo.seq_length = 120; % length of simulated traces in time steps


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSteps = 5.9;
trueParams.nSteps = nSteps;
trueParams.discrete_data_flag = 0;
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
mcmcInfo.PFResamplingFlag = 1;

mcmcInfo = setMCMCOptions(mcmcInfo);
if ~mcmcInfo.inferNStepsFlag
    mcmcInfo.nSteps = trueParams.nSteps;
end    

mcmcInfo = initializeInferenceArrays(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
mcmcInfo = inferenceWrapper(mcmcInfo);
toc

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ref_chain_ids = repelem(1:mcmcInfo.n_chains,mcmcInfo.n_temps_per_chain);
cmap = brewermap(mcmcInfo.n_chains,'Spectral');

close all
n_points = size(mcmcInfo.v_inf_array(:,2,:),1);
x_axis = mcmcInfo.update_increment*(1:n_points);

if mcmcInfo.nStates > 2
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
end
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

if mcmcInfo.inferNStepsFlag
    Mem_fig = figure;
    hold on
    for i = 1:mcmcInfo.n_chains_eff
        if mcmcInfo.refChainVec(i)
            plot(x_axis,mcmcInfo.n_steps_inf_array(:,i),'-o','Color',cmap(ref_chain_ids(i),:),'LineWidth',2);
        else
            plot(x_axis,mcmcInfo.n_steps_inf_array(:,i),'-','Color',[cmap(ref_chain_ids(i),:) .5],'LineWidth',1);
        end
    end
    plot(x_axis,repelem(mcmcInfo.trueParams.nSteps,n_points),'k--','LineWidth',2)

    ylabel('memory (n steps)')
    xlabel('MCMC steps')
end

% calculate first index to use
first_ind = mcmcInfo.burn_in/mcmcInfo.update_increment + 1;

if mcmcInfo.nStates>2
    v3_fig_hist = figure;
    v3_data = [];
    for i = find(mcmcInfo.refChainVec)
        v3_data = [v3_data mcmcInfo.v_inf_array(first_ind:end,3,i)'];        
    end
    histogram(v3_data);
    % plot(x_axis,repelem(mcmcInfo.trueParams.v(3),n_points),'k--','LineWidth',2)
    xlabel('state 3 emission rate')
    xlabel('probability')
end

v2_fig_hist = figure;
v2_data = [];
for i = find(mcmcInfo.refChainVec)
    v2_data = [v2_data mcmcInfo.v_inf_array(first_ind:end,2,i)'];        
end
histogram(v2_data);
% plot(x_axis,repelem(mcmcInfo.trueParams.v(3),n_points),'k--','LineWidth',2)
xlabel('state  emission rate')
ylabel('probability')

pon_fig_hist = figure;
hold on
a21_data = [];
for i = find(mcmcInfo.refChainVec)
    a21_data = [a21_data reshape(mcmcInfo.A_inf_array(2,1,first_ind:end,i),1,[])];        
end
h = histogram(a21_data,'Normalization','probability');
plot(repelem(mcmcInfo.trueParams.A(2,1),n_points),linspace(0,max(h.Values),n_points),'k--','LineWidth',2)
xlabel('pON')
ylabel('marginal posterior probability')

poff_fig_hist = figure;
hold on
a12_data = [];
for i = find(mcmcInfo.refChainVec)
    a12_data = [a12_data reshape(mcmcInfo.A_inf_array(1,2,first_ind:end,i),1,[])];        
end
h = histogram(a12_data,'Normalization','probability');
plot(repelem(mcmcInfo.trueParams.A(1,2),n_points),linspace(0,max(h.Values),n_points),'k--','LineWidth',2)
xlabel('pOff')
ylabel('marginal posterior probability')

%%
chain_id = 3;
A = mcmcInfo.A_curr(:,:,chain_id);
p0 = mcmcInfo.pi0_curr(chain_id,:);
v = mcmcInfo.v_curr(chain_id,:);
v(1) = 0;
pd_r = v*p0';

% get predicted mean fluo
f_mean_pd = sum(mcmcInfo.coeff_MS2(:,chain_id))*pd_r
% get actual mean fluo of trace predictions
f_mean_sim = mean(mean(mcmcInfo.sample_fluo(:,:,chain_id)))
% compare to mean fluo of actual traces
f_mean_true = mean(mcmcInfo.observed_fluo(:))