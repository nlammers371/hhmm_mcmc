% Script to validate ns MCMC method for burs parameter inference
clear 
close all force

addpath(genpath('../utilities'))

% make save diractor
DropboxFolder = 'S:\Nick\Dropbox (Personal)\';
if ~exist(DropboxFolder)
    DropboxFolder = 'C:\Users\nlamm\Dropbox (Personal)\';
end    
outPath = [DropboxFolder 'hhmm_MCMC_data\mcmc_v1_validation\'];
inf_files = dir([outPath '*.mat']);
master_struct = struct;
wb = waitbar(0,'loading MCMC results...');
for i = 1:length(inf_files)
    waitbar(i/length(inf_files),wb)
    load([outPath inf_files(i).name])
    % store key fields
    master_struct(i).sample_fluo = mcmcInfo.sample_fluo;
    master_struct(i).observed_fluo = mcmcInfo.observed_fluo;
    master_struct(i).trueParams = mcmcInfo.trueParams;
    master_struct(i).A_inf_array = mcmcInfo.A_inf_array;
    master_struct(i).v_inf_array = mcmcInfo.v_inf_array;
    master_struct(i).logL_array = mcmcInfo.logL_vec;
    master_struct(i).sigma_inf_array = mcmcInfo.sigma_inf_array;
    master_struct(i).n_traces = mcmcInfo.n_traces;
    master_struct(i).nSteps = mcmcInfo.nSteps;
    master_struct(i).nStates = mcmcInfo.nStates;
    clear mcmcInfo
end
delete(wb)

%% Organize results into groups
burn_in = 1e3;
stat_array = [[master_struct.nStates]' [master_struct.n_traces]' [master_struct.nSteps]'];
[stat_index,~,group_id_vec] = unique(stat_array,'rows');

for s = size(stat_index,1)-25
  
    % concatenate results from same group
    group_indices = find(group_id_vec==s);
    A_inf_array = [];
    v_inf_array = [];    
    logL_array = [];
    for g = 1:length(group_indices)
        A_inf_array = cat(4,A_inf_array,master_struct(group_indices(g)).A_inf_array(:,:,burn_in:end,:));
        v_inf_array = cat(3,v_inf_array,master_struct(group_indices(g)).v_inf_array(burn_in:end,:,:));
        logL_array = cat(2,logL_array,master_struct(group_indices(g)).logL_array(burn_in:end,:));
    end
    n_chains = size(v_inf_array,3);
    trueParams = master_struct(group_indices(1)).trueParams;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test different strategies for parameter estimation
    
    % (1) Take most likely chain
    logL_vec = mean(logL_array,1);
    [~ , si] = sort(logL_vec,'descend');
    
    [true_val_vec, param_mean_vec1, param_ste_vec1, param_names] = ...
                get_parameter_estimates(trueParams,A_inf_array(:,:,:,si(1)),v_inf_array(:,:,si(1)),0);
              
    % (2) Take average across all chains         
    [~, param_mean_vec2, param_ste_vec2, ~] = ...
                get_parameter_estimates(trueParams,A_inf_array,v_inf_array,0);
              
    % (3) Take average, removing outlier chains    
    [~, param_mean_vec3, param_ste_vec3, ~] = ...
                get_parameter_estimates(trueParams,A_inf_array,v_inf_array,1);
    
  
end  



