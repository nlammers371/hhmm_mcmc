% Script to validate ns MCMC method for burs parameter inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
if ~exist(DropboxFolder)
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
end    
runName = 'mcmc_validation_basic_mh';
% runName = 'run_mcmc_validation_v2_bootstrapping_';
outPath = [DropboxFolder 'hhmm_MCMC_data\' runName filesep];
figPath = ['../../fig/validation/' runName filesep];
mkdir(figPath)

inf_files = dir([outPath '*.mat']);
master_struct = struct;
wb = waitbar(0,'loading MCMC results...');
iter = 1;
for i = 1:length(inf_files)
    waitbar(i/length(inf_files),wb)
    try
        load([outPath inf_files(i).name])
        % store key fields
        master_struct(iter).sample_fluo = mcmcInfo.sample_fluo;
        master_struct(iter).observed_fluo = mcmcInfo.observed_fluo;
        master_struct(iter).trueParams = mcmcInfo.trueParams;
        master_struct(iter).A_inf_array = mcmcInfo.A_inf_array;
        master_struct(iter).err_flag_vec = mcmcInfo.err_flag_vec;
        master_struct(iter).v_inf_array = mcmcInfo.v_inf_array;
        master_struct(iter).logL_array = mcmcInfo.logL_vec;
        master_struct(iter).sigma_inf_array = mcmcInfo.sigma_inf_array;
        master_struct(iter).n_traces = mcmcInfo.n_traces_per_chain;
        master_struct(iter).nSteps = mcmcInfo.nSteps;
        master_struct(iter).nStates = mcmcInfo.nStates;
        clear mcmcInfo
        iter = iter + 1;
    catch
      % do nothing
    end
end
delete(wb)

% Organize results into groups
markerSize = 75;
burn_in = 2.5e3;
stat_array = [[master_struct.nStates]' [master_struct.n_traces]' [master_struct.nSteps]'];
[stat_index,~,group_id_vec] = unique(stat_array,'rows');
%%
wb = waitbar(0,'generating fit figures');
for s = 1:size(stat_index,1)
    waitbar(s/size(stat_index,1),wb)
    close all
    % concatenate results from same group
    group_indices = find(group_id_vec==s);
    n_traces = stat_index(s,2);
    nStates = stat_index(s,1);
    nSteps = stat_index(s,3);
    
    A_inf_array = [];
    v_inf_array = [];    
    logL_array = [];
    err_flag_vec = [];
    boot_flag_vec = [];
    for g = 1:length(group_indices)
        A_inf_array = cat(4,A_inf_array,master_struct(group_indices(g)).A_inf_array(:,:,burn_in:end,:));
        v_inf_array = cat(3,v_inf_array,master_struct(group_indices(g)).v_inf_array(burn_in:end,:,:));
        logL_array = cat(2,logL_array,master_struct(group_indices(g)).logL_array(burn_in:end,:));
        err_flag_vec = [err_flag_vec master_struct(group_indices(g)).err_flag_vec];
        sz = size(master_struct(group_indices(g)).A_inf_array(:,:,burn_in:end,:),4);
        boot_flag_vec = [boot_flag_vec repelem(g,sz)];
    end
    n_chains = size(v_inf_array,3);
    trueParams = master_struct(group_indices(1)).trueParams;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test different strategies for parameter estimation
    
    % (1) Take most likely chain
    logL_vec = mean(logL_array,1);
    [~ , si] = sort(logL_vec,'descend');
    
    [true_val_vec, param_mean_vec1, param_ste_vec1, param_names] = ...
                get_parameter_estimates(trueParams,A_inf_array(:,:,:,si(1)),v_inf_array(:,:,si(1)),0,false);
              
    % (2) Take average across all chains         
%     [~, param_mean_vec2, param_ste_vec2, ~] = ...
%                 get_parameter_estimates(trueParams,A_inf_array,v_inf_array,0,err_flag_vec);

    [~, param_mean_vec2, param_ste_vec2, ~] = ...
                get_parameter_estimates_boot(trueParams,A_inf_array,v_inf_array,nanmean(logL_array,1),boot_flag_vec,err_flag_vec);
              
    % (3) Take average, removing outlier chains    
    [~, param_mean_vec3, param_ste_vec3, ~] = ...
                get_parameter_estimates(trueParams,A_inf_array,v_inf_array,1,err_flag_vec);
              
              
    % generate basic results plot
    ind = 2;
    if trueParams.nStates == 3
        ind = 6;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Transition probability
    a_fig = figure('Visible','off','Position',[100 100 512 512]);
    cmap = brewermap([],'Set2');
    hold on
    xf1 = rand(1,ind)*0.1-0.05;
    xf2 = rand(1,ind)*0.1-0.05;
    xf3 = rand(1,ind)*0.1-0.05;
    
    % strat 1
    errorbar((1:ind)+xf1,param_mean_vec1(1:ind),param_ste_vec1(1:ind),'ok','Capsize',0)
    s1 = scatter((1:ind)+xf1,param_mean_vec1(1:ind),markerSize,'d','MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k');
    % strat 2
    errorbar((1:ind)+xf2,param_mean_vec2(1:ind),param_ste_vec2(1:ind),'ok','Capsize',0)
    s2 = scatter((1:ind)+xf2,param_mean_vec2(1:ind),markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');
    % strat 3
    errorbar((1:ind)+xf3,param_mean_vec3(1:ind),param_ste_vec3(1:ind),'ok','Capsize',0)
    s3 = scatter((1:ind)+xf3,param_mean_vec3(1:ind),markerSize,'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
    % ground truth
    s4 = scatter(1:ind,true_val_vec(1:ind),markerSize,'s','MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');
    
    grid on
    xlim([0.5 ind+1])
    set(gca,'xtick',1:ind,'xticklabel',param_names(1:ind))
    set(gca,'Fontsize',14)
    
    ylabel('transition probability')
    xlabel('parameter')
    legend([s1 s2 s3 s4], 'ML','ML (boot)','Mean','Truth','Location','northwest')
    
    title(['Transition prob results (K' num2str(nStates) ' W' num2str(nSteps) ' N' sprintf('%04d',100*n_traces) ')'])
    saveString = ['K' sprintf('%01d',nStates) '_W' sprintf('%03d',10*round(nSteps,1)) '_nt' sprintf('%03d',n_traces)];
    saveas(a_fig,[figPath 'a_fit_' saveString '.png'])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Emission rate
    dim = length(param_mean_vec1) - ind; 
    v_fig = figure('Visible','off','Position',[100 100 512 512]);
    cmap = brewermap([],'Set2');
    hold on
    xf1 = rand(1,dim)*0.1-0.05;
    xf2 = rand(1,dim)*0.1-0.05;
    xf3 = rand(1,dim)*0.1-0.05;
    
    % strat 1
    errorbar((1:dim)+xf1,param_mean_vec1(ind+1:end),param_ste_vec1(ind+1:end),'ok','Capsize',0)
    s1 = scatter((1:dim)+xf1,param_mean_vec1(ind+1:end),markerSize,'d','MarkerFaceColor',cmap(5,:),'MarkerEdgeColor','k');
    % strat 2
    errorbar((1:dim)+xf2,param_mean_vec2(ind+1:end),param_ste_vec2(ind+1:end),'ok','Capsize',0)
    s2 = scatter((1:dim)+xf2,param_mean_vec2(ind+1:end),markerSize,'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','k');
    % strat 3
    errorbar((1:dim)+xf3,param_mean_vec3(ind+1:end),param_ste_vec3(ind+1:end),'ok','Capsize',0)
    s3 = scatter((1:dim)+xf3,param_mean_vec3(ind+1:end),markerSize,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','k');
    % ground truth
    s4 = scatter(1:dim,true_val_vec(ind+1:end),markerSize,'s','MarkerFaceColor',cmap(8,:),'MarkerEdgeColor','k');
    
    grid on
    xlim([0.5 dim+1])
    set(gca,'xtick',1:dim,'xticklabel',param_names(ind+1:end))
    set(gca,'Fontsize',14)
    
    ylabel('emission rate (au)')
    xlabel('parameter')
    legend([s1 s2 s3 s4], 'ML','ML (boot)','Mean','Truth','Location','Southeast')
    
    title(['Emission rate results (K' num2str(nStates) ' W' num2str(nSteps) ' N' sprintf('%04d',100*n_traces) ')'])
    
    saveas(v_fig,[figPath 'v_fit_' saveString '.png'])
    
end 

delete(wb)