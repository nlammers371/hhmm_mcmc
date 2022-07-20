function mcmcInfo = resample_chains_PF(mcmcInfo)

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
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now conduct backward smoothing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sample_chains_bkd = NaN(size(sample_chains_fwd));
sample_fluo_bkd = NaN(size(sample_chains_fwd));
sample_fluo_prob_bkd = NaN(size(sample_chains_fwd,1)-1,n_chains,n_traces);
cp_state_array_bkd = NaN(size(cp_state_array_fwd));
cp_state_array_bkd(end,:,:) = cp_state_array_fwd(end,:,:);

sample_fluo_bkd(end,:,:) = sum(cp_init_tracker.*ms2_kernel,1); 
sample_chains_bkd(end-max_w+1:end,:,:) = sample_chains_fwd(end-max_w+1:end,:,:);    

% initialize tracker for state transitions
% tr_tracker_bkd = zeros(nStates, nStates, seq_len-1, n_traces);
% tr_tracker_bkd(:,:,end-max_w+2:end,:) = tr_tracker_fwd(:,:,end-max_w+2:end,:);

% add final fluo value
tic
% iterate through indices to sample
for i = seq_len-1:-1:1
    
    mcmcInfo.samp_index = i;        
%     post_time_index = i + 1;
    
    %%% resample based off of fluo from previous time point %%%
    
    % calculate fluorescence
    fluo_post = sum(cp_init_tracker.*ms2_kernel,1);           

    % calculate fluo probabilities 
    logL_fluo_full = -0.5*(((reshape(fluo_ref(i+1,:),1,1,[])-fluo_post)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
    probs_fluo = exp(logL_fluo_full-logsumexp(logL_fluo_full,2));

    sample_fluo_prob_bkd(i+1,:,:) = probs_fluo;
    
    % resample
    option_array = [zeros(1,1,n_traces) cumsum(probs_fluo,2)];
    rand_array = rand(1,n_chains,n_traces);    

    if i < seq_len-1
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
    
    if i >= max_w
        % extract resampled post states
        post_states = cp_state_tracker(end,:,:);

        % calculate linear indices 
        row_col_array_to = post_states + (mcmcInfo.row_ref-1)*nStates;% + mcmcInfo.chain_id_ref*nStates^2;

        % extract probabilities
        tr_probs = A(row_col_array_to);
        tr_probs = tr_probs ./ sum(tr_probs,1);

        % simulate jumps
        option_array = cumsum(tr_probs);
        rand_array_prop = repmat(rand(1,n_chains,n_traces),nStates,1);
        new_states = sum(rand_array_prop > option_array) + 1;
        sample_chains_bkd(i-max_w+1,:,:) = new_states;
                
%         tr_counts_temp = tr_tracker_bkd(:,:,i-max_w+1,:);
%         for n = 1:nStates
%             for m = 1:nStates
%                 tr_counts_temp(n,m,1,:) = reshape(sum(new_states==m & post_states==n,2),1,1,1,[]);
%             end
%         end
%         tr_tracker_bkd(:,:,i-max_w+1,:) = tr_counts_temp;
    
    else
        new_states = ones(size(new_states));
    end
                
    %%% update compound state tracker
    cp_state_tracker = circshift(cp_state_tracker,-1,1);
    cp_state_tracker(end,:,:) = new_states;        
    cp_state_array_bkd(i,:,:) = 1 + sum(cp_conv_array.*(cp_state_tracker-1),1);                   
    
    %%% initiation state tracker
    cp_init_tracker = circshift(cp_init_tracker,-1,1);
    if i >= max_w  
        v_lin_indices = 1 + n_chains*(new_states-1);
        cp_init_tracker(end,:,:) = mcmcInfo.v_curr(v_lin_indices);
    else               
        cp_init_tracker(end,:,:) = 0;
    end
    
    sample_fluo_bkd(i,:,:) = sum(cp_init_tracker.*ms2_kernel,1);   
end 
toc

%% calculate alpha (fwd) and beta (bkd) matrices (in naive space)
fwd_array = ones(seq_len,nStates,n_traces);
bkd_array = ones(seq_len,nStates,n_traces);
for i = 1:nStates 
    fwd_array(:,i,:) = mean(sample_chains_fwd==i,2);
    bkd_array(:,i,:) = mean(sample_chains_bkd==i,2);
end

% calculate joint state probability array <s>
joint_array = fwd_array.*bkd_array;
joint_array = joint_array ./ sum(joint_array,2);

%% 
ps_array = zeros(seq_len,nStates,n_traces);

n_vec = zeros(seq_len,n_traces);
tic
for t = 1:n_traces
    for i = 1:seq_len 
%         pss_slice = pss_array(:,:,i-1,t);        
        
        % get lists of unique compound states
        cpb = cp_state_array_bkd(i,:,t);
        [cpb_u,ia_b,ic_b] = unique(cpb);
        cpf = cp_state_array_fwd(i,:,t);
        [cpf_u,~,ic_f] = unique(cpf);
        
        % calculate probabilities for all compatible pairs        
        for s = 1:length(cpf_u)
            cp_from = cpf_u(s);            
            matching_states = ismember(cpb_u,cp_from);
            
            for m = find(matching_states)
              
                % get occurrence counts
                bc = sum(ic_b==m);
                fc = sum(ic_f==s);
                
                n_vec(i,t) = n_vec(i,t)+bc*fc;
            end
        end
%         n_vec(i-1,t) = sum(pss_slice(:));
%         if sum(pss_slice(:))~=0
%             pss_array(:,:,i-1,t) = pss_slice / sum(pss_slice(:));                
%         end
    end    
end

%% now calculate transition probabilities
pss_array = zeros(nStates,nStates,seq_len-1,n_traces);
% pss_array_norm = zeros(nStates,nStates,seq_len,n_traces);
n_vec = zeros(seq_len-1,n_traces);
tic
for t = 1:n_traces
    for i = 2:seq_len 
        pss_slice = pss_array(:,:,i-1,t);
        
        % get fluo probabilities
        pyi = sample_fluo_prob_bkd(i,:,t);
        
        % get lists of unique compound states
        cpb = cp_state_array_bkd(i,:,t);
        [cpb_u,ia_b,ic_b] = unique(cpb);
        cpf = cp_state_array_fwd(i-1,:,t);
        [cpf_u,~,ic_f] = unique(cpf);
        
        % calculate probabilities for all compatible pairs        
        for s = 1:length(cpf_u)
            cp_from = cpf_u(s);
            cp_to = allowed_from(cp_from, nStates, max_w);
            matching_states = ismember(cpb_u,cp_to);
            
            for m = find(matching_states)
                from_naive = digit(cp_from, 1, nStates, max_w);%mod(cp_from, nStates)+1;
                to_naive = digit(cpb_u(m), 1, nStates, max_w);%mod(cpb_u(m), nStates)+1;
                
                % get transition probability
                a_prob = A(to_naive,from_naive);
                
                % get occurrence counts
                bc = sum(ic_b==m);
                fc = sum(ic_f==s);
                
                pss_slice(to_naive,from_naive) = pss_slice(to_naive,from_naive) + a_prob*bc*fc*pyi(ia_b(m));
            end
        end
        n_vec(i-1,t) = sum(pss_slice(:));
        if sum(pss_slice(:))~=0
            pss_array(:,:,i-1,t) = pss_slice / sum(pss_slice(:));                
        end
    end    
end
toc
mcmcInfo.sample_chains = sample_chains_fwd;
mcmcInfo.pss_array = pss_array;
mcmcInfo.ps_array = joint_array;

em_time = toc;
if mcmcInfo.em_timer_flag
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end