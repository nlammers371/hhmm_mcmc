function mcmcInfo = resample_chains_PF_v2(mcmcInfo)

% This script resamples the microscopic promoter state for each extant
% chain using particle filtering methods that mimic classical
% forward-backward algorithm

% extract parameters
A = mcmcInfo.A_curr(:,:,1);
pi0 = mcmcInfo.pi0_curr(1,:);
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

% initial_state_vec = NaN(1,n_chains,n_traces);
% for n = 1:n_chains
%     initial_state_vec(1,n,:) = reshape(randsample(1:nStates,n_traces,true,pi0(n,:)),1,1,n_traces);    
% end 

% generate fluorescent noise lookup table 
sigma_ref = repmat(mcmcInfo.sigma_curr(1),1,n_chains,n_traces);
fluo_ref = permute(mcmcInfo.observed_fluo,[1,3,2]);

% generate lookup table for fluorescence values corresponding to each state
% at each position along the gene body
max_w = find(mcmcInfo.coeff_MS2(:,1)>0,1,'last');
ms2_kernel = mcmcInfo.coeff_MS2(1:max_w,:);

% generate array for conversion from longform to short CP state notation
cp_vec = (nStates).^(max_w-1:-1:0);
cp_conv_array = repmat(cp_vec',1,n_chains,n_traces);

% initialize array to keep track of compound states
cp_state_tracker = ones(max_w,n_chains,n_traces);
cp_init_tracker = zeros(max_w,n_chains,n_traces);
cp_state_array_fwd = zeros(seq_len,n_chains,n_traces);        
F_array_fwd_long = zeros(seq_len,n_chains,n_traces,nStates);        

% array to store cp state fluo likelihoods
sample_fluo_prob_fwd = zeros(seq_len,n_chains,n_traces);        

% initialize new "predicted fluo" array
sample_fluo_fwd = NaN(size(sample_chains_fwd));

% initialize tracker for state transitions
% tr_tracker_fwd = zeros(nStates, nStates, seq_len-1, n_traces);

tic
% iterate through indices to sample
for i = 1:seq_len
    
    mcmcInfo.samp_index = i;        
    prev_time_index = i - 1;
    
    %%% simulate next state transition %%%
    
    % calculate linear indices 
    if i > 1
        prev_state_array = sample_chains_fwd(prev_time_index,:,:);
    
        row_col_array_from = ...
                             (prev_state_array-1)*nStates + mcmcInfo.row_ref;    % mcmcInfo.chain_id_ref*nStates^2 + ...
    
        % extract probabilities
        tr_probs = A(row_col_array_from);
    
        % simulate jumps
        option_array = cumsum(tr_probs);
        rand_array_prop = repmat(rand(1,n_chains,n_traces),nStates,1);
        new_states = sum(rand_array_prop > option_array) + 1;
                        
    else                
        new_states = reshape(randsample(1:nStates,n_traces*n_chains,true,pi0),1,n_chains,n_traces);            
    end
    
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
    
    for t = 1:n_traces
      
        % resample
        rs_ids = discretize(rand_array(1,:,t),option_array(1,:,t));
        
        % resample arrays
        sample_chains_fwd(i,:,t) = sample_chains_fwd(i,rs_ids,t);
        sample_fluo_fwd(i,:,t) = fluo_next(1,rs_ids,t);
        if i > 1
            prev_state_array(1,:,t) = prev_state_array(1,rs_ids,t);
        end
        
        % resample compound state trackers
        cp_state_tracker(:,:,t) = cp_state_tracker(:,rs_ids,t);
        cp_init_tracker(:,:,t) = cp_init_tracker(:,rs_ids,t);
        
        % add fluo probabilities
        sample_fluo_prob_fwd(i,:,t) = probs_fluo(1, rs_ids, t);
    end
    
    % convert cp states to integer representation   
    cp_state_array_fwd(i,:,:) = 1 + sum(cp_conv_array.*(cp_state_tracker-1),1);    
    for k = 1:nStates
%         test(k,:,:) = sum(ms2_kernel.*(cp_state_tracker==k),1);
        F_array_fwd_long(i,:,:,k) = sum(ms2_kernel.*(cp_state_tracker==k),1);
    end
%     if i > 1
%         tr_counts_temp = tr_tracker_fwd(:,:,i-1,:);
%         for n = 1:nStates
%             for m = 1:nStates
%                 tr_counts_temp(n,m,1,:) = reshape(sum(prev_state_array==m & sample_chains_fwd(i,:,t)==n,2),1,1,1,[]);
%             end
%         end
%         tr_tracker_fwd(:,:,i-1,:) = tr_counts_temp;
%     end
end 
% toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now conduct backward smoothing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample_chains_bkd = NaN(size(sample_chains_fwd));
% sample_fluo_bkd = NaN(size(sample_chains_fwd));
% sample_fluo_prob_bkd = NaN(size(sample_chains_fwd,1)-1,n_chains,n_traces);
% cp_state_array_bkd = NaN(size(cp_state_array_fwd));
% cp_state_array_bkd(end,:,:) = cp_state_array_fwd(end,:,:);
bkd_prior = 1/n_chains^2;
sample_logL_bkd = ones(size(cp_state_array_fwd)) * log(bkd_prior);
sample_logL_bkd(end,:,:) = log(1/n_chains); 

% initialize tracker for state transitions

% tic
% iterate through indices to sample
for i = seq_len-1:-1:1
    
    mcmcInfo.samp_index = i;        
    
    %%% resample based off of fluo from previous time point %%%
    
    % calculate fluorescence
    fluo_post = sum(cp_init_tracker.*ms2_kernel,1);           

    % calculate fluo probabilities 
    logL_fluo_full = -0.5*(((reshape(fluo_ref(i+1,:),1,1,[])-fluo_post)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));       
    fluo_prob = exp(logL_fluo_full);
    
    % get previous prob
    logL_prev = sample_logL_bkd(i+1,:,:);
%     logL_prev = logL_prev - logsumexp(logL_prev,2);
    prob_prev = exp(logL_prev);
    
    % resample
    cp_to_array = cp_state_array_fwd(i+1,:,:);
    cp_from_array = cp_state_array_fwd(i,:,:);    
    
    for t = 1:n_traces                
        
        logL_slice = sample_logL_bkd(i,:,t);
        
        cpb = cp_to_array(1,:,t);
        [cpb_u,ia_b,~] = unique(cpb);
        
        cpf = cp_from_array(1,:,t);
        [cpf_u,~,ic_f] = unique(cpf);
        
        % get array of possible states that could transition into i+1
        % states
        from_array = allowed_to_array(cpb_u, nStates, max_w);
        
        % generate from array of matching size
        to_array = repmat(cpb_u,nStates,1);
        
        % get first digits of each
        from_n = digit(from_array, 1, nStates, max_w);
        to_n = digit(to_array, 1, nStates, max_w);
        
        % calculate transition probabilities
        a_lin_indices = to_n + (from_n-1)*nStates; 
        tr_probs = A(a_lin_indices);
        
        % get total probability
        bkd_probs = tr_probs.*fluo_prob(ia_b').*prob_prev(ia_b');
        
        % now assign to prev states        
        for s = 1:length(cpf_u)
            cp_from = cpf_u(s);
            if any(cp_from==from_array(:))
                logL_slice(1,ic_f==s) = log(sum(bkd_probs(cp_from==from_array)));
            end
        end                        
        
        % resample
        sample_logL_bkd(i,:,t) = logL_slice - logsumexp(logL_slice,2);
        if any(isnan(sample_logL_bkd(i,:,t)))
            error('wtf')
        end
    end
    
end 
% toc

mcmcInfo.sample_chains = sample_chains_fwd;
mcmcInfo.sample_logL_bkd = sample_logL_bkd;
mcmcInfo.sample_fluo_prob = sample_fluo_prob_fwd;
mcmcInfo.sample_fluo = sample_fluo_fwd;
mcmcInfo.cp_state_array = cp_state_array_fwd;

%% Calculate transition probability array

pss_array = zeros(nStates,nStates,seq_len-1,n_traces);
% tic
for t = 1:n_traces
    for i = 2:seq_len
        pss_slice = pss_array(:,:,i-1,t);
        
        % extract relevant array slices
        cp_from_array = cp_state_array_fwd(i-1,:,t); 
        cp_to_array = cp_state_array_fwd(i,:,t); 
        bkd_prob_array = exp(sample_logL_bkd(i,:,t));
        fluo_prob_array = sample_fluo_prob_fwd(i,:,t);
        
        % get array of all states that could transition to s_i
%         [cpf_u,ia_f,~] = unique(cp_from_array);
        [cpb_u,ia_b,~] = unique(cp_to_array);
        from_array = allowed_to_array(cpb_u, nStates, max_w);
        
        % get first digits of each
        from_n = digit(cpb_u, 2, nStates, max_w);
        to_n = digit(cpb_u, 1, nStates, max_w);
        
        % calculate transition probabilities
        a_lin_indices = to_n + (from_n-1)*nStates; 
        tr_probs = A(a_lin_indices);
        
        % get count of s_{i-1} states transitioning to each s_i
        n_si1_vec = zeros(1,length(cpb_u));
        for s = 1:length(cpb_u)
            n_si1_vec(s) = sum(ismember(cp_from_array,from_array(:,s)));
        end
        
        % put it all together
        total_probs = fluo_prob_array(ia_b).*tr_probs.*n_si1_vec; %bkd_prob_array(ia_b).*
        for k = 1:nStates
            for l = 1:nStates
                pss_slice(l,k) = sum(total_probs(to_n == l & from_n == k));
            end
        end
        pss_array(:,:,i-1,t) = pss_slice / sum(pss_slice(:));
    end       
end 


%% Now calculate p_s matrix

ps_array_full = exp(sample_logL_bkd).*F_array_fwd_long;
F_array = permute(sum(ps_array_full,2)./sum(exp(sample_logL_bkd),2),[1 4 3 2]);
if any(isnan(F_array(:)))
    error(';wtf')
end    
toc

% record
mcmcInfo.pss_array = pss_array;
mcmcInfo.F_array = F_array;

%%
em_time = toc;
if mcmcInfo.em_timer_flag
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end