% Script to make illustrative inference figures for 2 state case
clear 
close all

addpath(genpath('../utilities'))

% make path to write figures
figPath = '../../fig/em_run_time/';
mkdir(figPath)

% initialize info structure
trueParams = setParamsBasic3state;
% trueParams.sigma = 1;

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
rng(458)
%%
n_mcmc_steps = 3e2;
mcmcInfo.n_mcmc_steps = n_mcmc_steps; % number of MCMC steps (need to add convergence criteria)
mcmcInfo.burn_in = 100;
n_chains = 1;
mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run
mcmcInfo.n_reps = 2; % number of chain state resampling passes per inference step

% characteristics of simulated data
mcmcInfo.n_traces = 1;
mcmcInfo.seq_length = 240; % length of simulated traces in time steps
mcmcInfo.em_timer_flag = true;
mcmcInfo.NumWorkers = 4;
mcmcInfo.inferNStepsFlag = false;
mcmcInfo.trueParams = trueParams;
mem_vec = 3:15;

% initialize array
em_time_vec = NaN(size(mem_vec));
master_struct = struct;
%%
initializePool(mcmcInfo,0);

for m = 1:length(mem_vec)
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set MCMC options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nSteps = mem_vec(m);
    trueParams = mcmcInfo.trueParams;
    trueParams.nSteps = nSteps;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simulate data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfoTemp = mcmcInfo;
    trueParams.n_traces = mcmcInfoTemp.n_traces;
    trueParams.seq_length = mcmcInfoTemp.seq_length;
    trueParams = generateSimulatedData(trueParams);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize inference arrays and variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfoTemp.nStates = trueParams.nStates;    
    mcmcInfoTemp.alpha_frac = trueParams.alpha_frac;
    mcmcInfoTemp.observed_fluo = trueParams.observed_fluo;

    mcmcInfoTemp = setMCMCOptions(mcmcInfoTemp);
    mcmcInfoTemp.nSteps = trueParams.nSteps;
    mcmcInfoTemp.nStepsMax = ceil(trueParams.nSteps)+1;
    
    mcmcInfoTemp = initializeInferenceArrays(mcmcInfoTemp);
    mcmcInfoTemp = initializeVariablesBasicRandom(mcmcInfoTemp);

    mcmcInfoTemp.save_trace_results = true;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % conduct inference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tic
    mcmcInfoTemp = inferenceWrapper(mcmcInfoTemp);
%     toc
    master_struct(m).mcmcInfo = mcmcInfoTemp;
    master_struct(m).em_time_mcmc_vec = mcmcInfoTemp.em_time_vec;
    em_time_vec(m) = nanmean(mcmcInfoTemp.em_time_vec);
end

%% cpHMM
% call local EM function;
fluo = {master_struct(1).mcmcInfo.observed_fluo(:,1)'};
nStates = 3;
A_log = log(trueParams.A);
pi0_log = log(trueParams.pi0');
kappa = 1;%trueParams.alpha;
n_steps_max_vec = floor(logspace(2,0,length(mem_vec)));
eps = 1e-50;
noise = trueParams.sigma;
v = trueParams.v';

em_time_vec_cpHMM = NaN(size(mem_vec));

% should take ~20 minutes total to run
for m = 1:find(mem_vec<10,1,'last') % NL: getting an error past 9...
    n_steps_max = n_steps_max_vec(m);
    tic
    nSteps = mem_vec(m);
    local_em_outputs = local_em_MS2_reduced_memory_timer (fluo, ...
                  v, noise, pi0_log, A_log, nStates, nSteps, kappa, n_steps_max, eps);
    master_struct(m).local_em = local_em_outputs;
    master_struct(m).em_time_vec_cpHMM = local_em_outputs.em_time_vec;
    em_time_vec_cpHMM(m) = nanmean(local_em_outputs.em_time_vec);
    toc
end                
            
%% Make figures 

% perform simple fit to extrapoltate cpHMM performance
last_ind = find(~isnan(em_time_vec_cpHMM),1,'last');
time_fun = @(x) x(1)*x(2).^mem_vec(2:last_ind);
loss_fun = @(x) log10(time_fun(x)) - log10(em_time_vec_cpHMM(2:last_ind));
x_fit = lsqnonlin(loss_fun,[0.01 3],[0 3],[10 3]);
cpHMM_time_extrap = x_fit(1)*x_fit(2).^mem_vec;

close all
bkg_color = [228,221,209]/255;

% calculate std
em_time_std = NaN(size(mem_vec));
em_time_std_cpHMM = NaN(size(mem_vec));
for m = 1:length(mem_vec)
    em_time_std(m) = nanstd(master_struct(m).em_time_mcmc_vec);
    em_time_std_cpHMM(m) = nanstd(master_struct(m).em_time_vec_cpHMM);
end    

em_fig = figure;
cmap = brewermap([],'Set2');
hold on

p = plot(mem_vec,cpHMM_time_extrap,'-.','Color',cmap(4,:),'LineWidth',2);
errorbar(mem_vec,em_time_vec,em_time_std,'Color','k','CapSize',0)
errorbar(mem_vec,em_time_vec_cpHMM,em_time_std_cpHMM,'Color','k','CapSize',0)

s(1) = scatter(mem_vec,em_time_vec,'MarkerEdgeColor','k','MarkerFaceColor',cmap(5,:));
s(2) = scatter(mem_vec,em_time_vec_cpHMM,'MarkerEdgeColor','k','MarkerFaceColor',cmap(4,:));

xlabel('memory (N time steps)')
ylabel('EM step time (seconds)')
xlim([3 15])
set(gca,'Fontsize',14)
% ylim([-0.1 4.5])
legend([s p],'MCMC','cpHMM','extrapolation','Location','northwest')
set(gca,'Color',bkg_color) 
box on
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.XAxis(1).Color = 'k';
set(gca,'yscale','log')
em_fig.InvertHardcopy = 'off';
set(gcf,'color','w');

saveas(em_fig,[figPath 'em_time_fig.png'])
saveas(em_fig,[figPath 'em_time_fig.pdf'])
