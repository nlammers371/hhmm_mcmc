function mcmcInfo = resample_chains_PF(mcmcInfo)

% This script resamples the microscopic promoter state for each extant
% chain using particle filtering methods that mimic classical
% forward-backward algorithm

% extract parameters
A = mcmcInfo.A_curr;
pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains_eff;
nStepsMax = mcmcInfo.nStepsMax;
seq_len = mcmcInfo.seq_length;

% generate reference vector
mcmcInfo.chain_id_ref = 0:n_chains-1;
mcmcInfo.trace_id_ref = reshape(0:n_traces-1,1,1,[]);
mcmcInfo.row_ref = (1:nStates)';
mcmcInfo.step_ref = (-nStepsMax+1:nStepsMax-1)';

% generate temporary chain array that includes post and prior states to
% make resampling easier below
sample_chains_fwd = NaN(size(mcmcInfo.sample_chains));
% sample_emissions_temp = NaN(size(mcmcInfo.sample_chains));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct forward particle simulations first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_state_vec = NaN(1,n_chains,n_traces);
for n = 1:n_chains
    initial_state_vec(1,n,:) = reshape(randsample(1:nStates,n_traces,true,pi0(n,:)),1,1,n_traces);    
end 

% generate fluorescent noise lookup table 
sigma_ref = repmat(mcmcInfo.sigma_curr',1,1,n_traces);
fluo_ref = permute(mcmcInfo.observed_fluo,[1,3,2]);

% generate lookup table for fluorescence values corresponding to each state
% at each position along the gene body
max_w = find(mcmcInfo.coeff_MS2(:,1)>0,1,'last');
ms2_kernel = mcmcInfo.coeff_MS2(1:max_w,:);

% generate ref array for gene position
gp_ref_array = repmat((1:max_w)',1,n_chains,n_traces);

% initialize array to keep track of compound states
cp_state_tracker = zeros(max_w,n_chains,n_traces);
cp_init_tracker = zeros(max_w,n_chains,n_traces);
        
% initialize new "predicted fluo" array
fluo_pd_array_fwd = NaN(size(sample_chains_fwd));

tic
% iterate through indices to sample
for i = 1:seq_len
    
    mcmcInfo.samp_index = i;        
    prev_time_index = i - 1;
    
    %%% simulate next state transition %%%
    
    % calculate linear indices 
    if i > 1
        prev_state_array = sample_chains_fwd(prev_time_index,:,:);
    else
        prev_state_array = initial_state_vec;
    end
    row_col_array_from = mcmcInfo.chain_id_ref*nStates^2 + (prev_state_array-1)*nStates+mcmcInfo.row_ref;    
    
    % extract probabilities
    tr_probs = A(row_col_array_from);
    
    % simulate jumps
    option_array = cumsum(tr_probs);
    rand_array_prop = repmat(rand(1,n_chains,n_traces),nStates,1);
    new_states = sum(rand_array_prop > option_array) + 1;
    sample_chains_fwd(i,:,:) = new_states;    
    
    %%% update compound state tracker
    cp_state_tracker = circshift(cp_state_tracker,1,1);
    cp_state_tracker(1,:,:) = new_states;        
    
    %%% initiation state tracker
    cp_init_tracker = circshift(cp_init_tracker,1,1);
    v_lin_indices = (mcmcInfo.chain_id_ref+1) + n_chains*(new_states-1);
    cp_init_tracker(1,:,:) = mcmcInfo.v_curr(v_lin_indices);
    
    %%% calculate fluorescence
    fluo_next = sum(cp_init_tracker.*ms2_kernel,1);    
    
    %%% resample based on difference from observed fluo %%%
    
    % calculate fluo probabilities 
    logL_fluo_full = -0.5*(((reshape(fluo_ref(i,:),1,1,[])-fluo_next)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
    probs_fluo = exp(logL_fluo_full-logsumexp(logL_fluo_full,2));    
       
    % resample
    option_array = [zeros(1,1,n_traces) cumsum(probs_fluo,2)];
    rand_array = rand(1,n_chains,n_traces);    
%     rs_ids = NaN(size(probs_fluo));
    for t = 1:n_traces
        % resample
        rs_ids = discretize(rand_array(1,:,t),option_array(1,:,t));
        
        % resample arrays
        sample_chains_fwd(i,:,t) = sample_chains_fwd(i,rs_ids,t);
        fluo_pd_array_fwd(i,:,t) = fluo_next(1,rs_ids,t);
        
        % resample compound state trackers
        cp_state_tracker(:,:,t) = cp_state_tracker(:,rs_ids,t);
        cp_init_tracker(:,:,t) = cp_init_tracker(:,rs_ids,t);
    end
                    
end 
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now conduct backward smoothing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_chains_bkd = NaN(size(sample_chains_fwd));
sample_fluo_bkd = NaN(size(sample_chains_fwd));

% initialize last max_w falues with results from forward sweep

tic
% iterate through indices to sample
for i = seq_len-max_w:-1:1
    
    mcmcInfo.samp_index = i;        
%     post_time_index = i + 1;
    
    %%% resample based off of fluo from previous time point %%%
    
    % calculate fluorescence
    fluo_post = sum(cp_init_tracker.*ms2_kernel,1); 
    sample_fluo_bkd(i+max_w,:,:) = fluo_post;
    
    if i < seq_len-max_w        

        % calculate fluo probabilities 
        logL_fluo_full = -0.5*(((reshape(fluo_ref(i+max_w,:),1,1,[])-fluo_post)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
        probs_fluo = exp(logL_fluo_full-logsumexp(logL_fluo_full,2));

        % resample
        option_array = [zeros(1,1,n_traces) cumsum(probs_fluo,2)];
        rand_array = rand(1,n_chains,n_traces);    

        for t = 1:n_traces
            % resample
            rs_ids = discretize(rand_array(1,:,t),option_array(1,:,t));

            % resample compound state trackers
            cp_state_tracker(:,:,t) = cp_state_tracker(:,rs_ids,t);
            cp_init_tracker(:,:,t) = cp_init_tracker(:,rs_ids,t);                

        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% simulate (reverse) jumps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract resampled post states
    post_states = cp_state_tracker(end,:,:);
    
    % calculate linear indices 
    row_col_array_to = mcmcInfo.chain_id_ref*nStates^2 + (mcmcInfo.row_ref-1)*nStates + post_states;
    
    % extract probabilities
    tr_probs = A(row_col_array_to);
    tr_probs = tr_probs ./ sum(tr_probs,1);
    
    % simulate jumps
    option_array = cumsum(tr_probs);
    rand_array_prop = repmat(rand(1,n_chains,n_traces),nStates,1);
    new_states = sum(rand_array_prop > option_array) + 1;
    sample_chains_bkd(i,:,:) = new_states;
        
    %%% update compound state tracker
    cp_state_tracker = circshift(cp_state_tracker,-1,1);
    cp_state_tracker(end,:,:) = new_states;        
    
    %%% initiation state tracker
    cp_init_tracker = circshift(cp_init_tracker,-1,1);
    v_lin_indices = (mcmcInfo.chain_id_ref+1) + n_chains*(new_states-1);
    cp_init_tracker(end,:,:) = mcmcInfo.v_curr(v_lin_indices);
                    
end 
toc

%% calculate a and b matrices
fwd_array = ones(seq_len,nStates,n_traces);
bkd_array = ones(seq_len,nStates,n_traces);
for i = 1:nStates 
    fwd_array(:,i,:) = mean(sample_chains_fwd==i,2);
    bkd_array(:,i,:) = mean(sample_chains_bkd==i,2);
end
bkd_array(end-max_w+1:end,:,:) = 0.5;

joint_array = fwd_array.*bkd_array;
joint_array = joint_array ./ sum(joint_array,2);

mcmcInfo.sample_chains = sample_chains_fwd(2:end-1,:,:);
em_time = toc;
if mcmcInfo.em_timer_flag
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end