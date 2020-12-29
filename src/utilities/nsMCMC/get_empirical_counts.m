function mcmcInfo = get_empirical_counts(mcmcInfo)
    
    % extract parameters
    A_curr = mcmcInfo.A_curr;
    A_log_full = log(mcmcInfo.A_curr);
    pi0_log_full = log(mcmcInfo.pi0_curr);
    nStates = size(A_curr,1);
    nSteps = mcmcInfo.nSteps;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    seq_length = mcmcInfo.seq_length;
    
    % initialize count arrays 
    mcmcInfo.transition_count_array = zeros(nStates,nStates,n_chains);
    mcmcInfo.state_count_array = zeros(n_chains,nStates);
    mcmcInfo.pi0_count_array = zeros(n_chains,nStates);
    
    % initialize likelihood array
    mcmcInfo.trace_logL_array = -Inf(seq_length,n_chains,n_traces);
    mcmcInfo.trace_logL_vec = -Inf(1,n_chains,n_traces);
    
    % generate kernel for fluo logL attribution
    fluo_logL_kernel = ones(nSteps,1)/nSteps;
    
    % iterate through traces
    for n = 1:n_chains
        
        if mcmcInfo.par_chain_flag
            A_log = A_log_full(:,:,n);
            pi0_log = pi0_log_full(n,:);
        else
            A_log = A_log_full;
            pi0_log = pi0_log_full;
        end
        
        sample_chains_slice = permute(mcmcInfo.sample_chains(:,n,:),[1 3 2]);
        logL_transition = -Inf(size(sample_chains_slice));
        
        % calculate initial state PDF        
        [GC,GR] = groupcounts(sample_chains_slice(1,:)');
        mcmcInfo.pi0_count_array(n,GR) = ...
                          mcmcInfo.pi0_count_array(n,GR) + GC';
                        
        % get likelihood
        logL_transition(1,:) = pi0_log(sample_chains_slice(1,:));        
        
        % calculate transition and occupancy statistics        
        for t = 1:seq_length-1
            % get transition counts
            from = sample_chains_slice(t,:);
            to = sample_chains_slice(t+1,:);
            lin_indices_tr = sub2ind(size(mcmcInfo.transition_count_array),...
                                      to,from,repelem(n,size(sample_chains_slice,2)));
            [GC,GR] = groupcounts(lin_indices_tr');
            mcmcInfo.transition_count_array(GR) = ...
                                  mcmcInfo.transition_count_array(GR) + GC;
                                
            % get likelihood
            logL_transition(t+1,:) = A_log(lin_indices_tr-(n-1)*nStates^2);
                                    
        end
               
        % get state occupancy counts             
        [GC,GR] = groupcounts(sample_chains_slice(:));
        mcmcInfo.state_count_array(n,GR) = GC;        
        
        % calculate fluo probabilities
        logL_fluo = -.5*((mcmcInfo.observed_fluo - permute(mcmcInfo.sample_fluo(:,n,:),...
              [1 3 2]))./mcmcInfo.sigma_curr(n)).^2 - log(sqrt(2*pi)*mcmcInfo.sigma_curr(n));
        logL_fluo_conv = convn(fluo_logL_kernel,logL_fluo,'full');  
        logL_fluo_conv = logL_fluo_conv(1:end-nSteps+1,:,:);
        
        % combine probability components
        mcmcInfo.trace_logL_array(:,n,:) = logL_fluo_conv + logL_transition;
        
        % get total trace likelihoods
        mcmcInfo.trace_logL_vec(1,n,:) = mean(mcmcInfo.trace_logL_array(:,n,:));
    end
    
    mcmcInfo.logL_vec(mcmcInfo.step) = mean(mcmcInfo.trace_logL_vec(:));
    
    % resample based on chain likelihood
    chain_weights = exp(mean(mcmcInfo.trace_logL_vec,3)/10); %NL: factor of 10 is an add hoc temperature variable
    chain_indices = randsample(1:n_chains,n_chains,true,chain_weights);
    
    mcmcInfo.sample_chains = mcmcInfo.sample_chains(:,chain_indices,:);
    mcmcInfo.transition_count_array = mcmcInfo.transition_count_array(:,:,chain_indices);