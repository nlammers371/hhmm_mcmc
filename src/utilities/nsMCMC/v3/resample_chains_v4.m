function mcmcInfo = resample_chains_v4(mcmcInfo)

% This script resamples the microscopic promoter state for each extant
% chain in an asynchronous manner. Linear indexing is used so that we can
% randomize the order in which each chain is resampled wrspt its neighbors.
% This should hopefully improve mixing time

% extract parameters
A_log = log(mcmcInfo.A_curr);
pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains;
n_temps = mcmcInfo.n_temps_per_chain;
nSteps = mcmcInfo.nSteps;

% generate random sampling orders
n_reps = mcmcInfo.n_reps;
sample_indices = cell(1,n_traces);
for m = 1:n_traces
    seq_length = length(mcmcInfo.observed_fluo{m});
    n_samples = seq_length*n_chains*n_temps;
%     for n = 1:n_chains
%         for t = 1:n_temps
%             temp_array(:,n,t) = randsample(repelem(1:seq_length,n_reps),n_reps*seq_length,false);
%         end
%     end
    sample_indices{m} = reshape(randsample(1:seq_length,n_samples,true),seq_length,n_chains,n_temps);%temp_array;
end

% generate temporary chain array that includes post and prior states to
% make resampling easier below
sample_chains_tr_probs = cell(1,n_traces);

for m = 1:n_traces
    seq_length = length(mcmcInfo.observed_fluo{m});
    % additional helper variables
    seq_len_temp(m) = seq_length+2;
    seq_len_dummy(m) = seq_length+2*(nSteps-1);
    % generate temp array
    sample_chains_tr_probs{m} = NaN(seq_length+2,n_chains,n_temps);    
    sample_chains_tr_probs{m}(2:end-1,:,:) = mcmcInfo.sample_chains{m};
end

% initialize post and prior states randomly using pi0 PDF
for m = 1:n_traces
    for n = 1:n_chains
        sample_chains_tr_probs{m}(1,n,:) = reshape(randsample(1:nStates,n_temps,true,pi0(n,:)),1,1,n_temps);
        sample_chains_tr_probs{m}(end,n,:) = reshape(randsample(1:nStates,n_temps,true,pi0(n,:)),1,1,n_temps);
    end     
end

% second temporary array used for fluorescence sampling
sample_chains_fluo = cell(1,n_temps);
observed_fluo_fluo = cell(1,n_temps);
% observed_fluo_dummy2 = cell(1,n_temps);
% sample_fluo_dummy = cell(1,n_temps);
for m = 1:n_traces
    sample_chains_fluo{m} = cat(1,ones(nSteps-1,n_chains,n_temps),mcmcInfo.sample_chains{m},ones(nSteps-1,n_chains,n_temps));
    observed_fluo_fluo{m} = cat(1,zeros(nSteps-1,1),mcmcInfo.observed_fluo{m},zeros(nSteps-1,1));

%     observed_fluo_dummy2{m} = cat(1,mcmcInfo.observed_fluo{m},zeros(nSteps-1,1)); % NL: is this used?
%     sample_fluo_dummy{m} = cat(1,mcmcInfo.sample_fluo{m},zeros(nSteps-1,n_chains,n_temps));
end


% Extract helper vectors
chain_id_ref = mcmcInfo.chain_id_ref;
temp_id_ref = mcmcInfo.temp_id_ref;
a_row_ref = mcmcInfo.row_ref;
temperingFlag = mcmcInfo.temperingFlag; 
tempGradVec = reshape(mcmcInfo.tempGradVec,1,1,[]);

% extract sample chains themselves
sample_chains_temp = mcmcInfo.sample_chains;
tic;
% iterate through indices to sample
for m = 1:n_traces
  
    % preallocate helper array for indexing
    index_helper1 = chain_id_ref*seq_len_temp(m) + temp_id_ref*seq_len_temp(m)*n_chains;
    index_helper2 = chain_id_ref*nStates^2 + (a_row_ref-1)*nStates;
    index_helper3 = chain_id_ref*seq_len_dummy(m) + temp_id_ref*seq_len_dummy(m)*n_chains;
    
    for i = 1:seq_length*n_reps  

%         mcmcInfo.indArray = sample_indices{m}(i,:,:);    
        prev_time_index_array = sample_indices{m}(i,:,:);    
        prev_lin_index_array = prev_time_index_array + index_helper1;%mcmcInfo.chain_id_ref*seq_len_temp + mcmcInfo.trace_id_ref*seq_len_temp*n_chains;%ind_helper;
        post_lin_index_array = prev_lin_index_array + 2;

        %%% previous state %%%
        % calculate linear indices 
        prev_state_array = sample_chains_tr_probs{m}(prev_lin_index_array);
        row_col_array_from = chain_id_ref*nStates^2 + (prev_state_array-1)*nStates+a_row_ref;    

        % extract probabilities
        prev_probs_log = A_log(row_col_array_from);

        %%% following state %%%
        % calculate linear indices 
        post_state_array = sample_chains_tr_probs{m}(post_lin_index_array);  
        row_col_array_to = index_helper2 + post_state_array;

        % extract probabilities
        post_probs_log = A_log(row_col_array_to);

        % combine
        logL_tr = prev_probs_log + post_probs_log;

        %%% calculate fluorescence probability component
        % calculate fluo error term       

        logL_fluo = calculate_fluo_logL_v4(mcmcInfo,sample_chains_fluo{m},observed_fluo_fluo{m},prev_time_index_array);           

        %%% put everything together
        total_log_likelihoods = logL_tr + logL_fluo;
        total_log_likelihoods = total_log_likelihoods - logsumexp(total_log_likelihoods,1);
        
        % apply differential temperature correction if appropirate
        if temperingFlag 
            total_log_likelihoods = total_log_likelihoods ./ tempGradVec;
        end
        total_likelihoods = exp(total_log_likelihoods);    

        %%% draw new samples
        option_array = cumsum(total_likelihoods);
        rand_array = repmat(rand(1,n_chains,n_temps),nStates,1);
        sample_chains_tr_probs{m}(prev_lin_index_array+1) = sum(rand_array > option_array) + 1;

        % pass along to dummy array
        curr_lin_index_array = (prev_time_index_array+nSteps-1) + index_helper3;
        sample_chains_fluo{m}(curr_lin_index_array) = sample_chains_tr_probs{m}(prev_lin_index_array+1);

    end 
    sample_chains_temp{m} = sample_chains_tr_probs{m}(2:end-1,:,:);
end
toc;
mcmcInfo.sample_chains = sample_chains_temp;