function mcmcInfo = resample_chains(mcmcInfo)

% extract parameters
A = mcmcInfo.A_curr;
A_log = log(mcmcInfo.A_curr);

nStates = size(A,1);
n_traces = mcmcInfo.n_traces;
nSteps = mcmcInfo.nSteps;
n_chains = mcmcInfo.n_chains;
seq_length = mcmcInfo.seq_length;
state_options = 1:nStates;
state_ref = repmat(reshape(state_options,1,1,[]),1,n_chains);
coeff_MS2 = mcmcInfo.coeff_MS2;
sigma = mcmcInfo.noise;
sigma_log = log(sigma);

% initialize
sample_chains_curr = mcmcInfo.sample_chains;

% interate through traces
for n = 1:n_traces
    % extract true fluorescence
    true_fluo = mcmcInfo.masterSimStruct(n).fluo_MS2';
    
    % slice sampling array
    sample_chains_slice = sample_chains_curr(:,:,n);
    
    % determine order of samplling
    sample_indices = randsample(1:seq_length,seq_length,false);
    
    % iterate through indices to sample
    for ind = sample_indices
        %%% calculate forward-backward state probabilities
        prev_probs_log = 0;                  
        if ind > 1
            prev_probs_log = A_log(:,sample_chains_slice(ind-1,:));            
        end
        post_probs_log = 0;        
        if ind < seq_length
            post_probs_log = A_log(sample_chains_slice(ind+1,:),:);            
        end
        % combine
        tr_probs_log = prev_probs_log + post_probs_log';
        
        %%% calculate fluorescence probability component
        
        % get predicted fluo fragment
        postInd = min([seq_length,ind+nSteps-1]);
        prevInd = max([1,ind-nSteps]);
        chain_fragment = repmat(sample_chains_slice(prevInd:postInd,:),1,1,nStates);        
        chain_fragment = vertcat(chain_fragment(1:ind-prevInd,:,:),state_ref,...
                                      chain_fragment(ind-prevInd+1:postInd-prevInd,:,:));
        fluo_fragment = convn(chain_fragment,coeff_MS2,'full');        
        fluo_fragment = fluo_fragment(ind-prevInd+1:end-nSteps+1,:,:);
        
        % calculate fluo error term
        log_fluo_diffs = reshape(-sum(((true_fluo(ind:postInd)-fluo_fragment)-sigma_log).^2,1),nStates,[]);
        
        %%% put everything together
        total_log_likelihoods = log_fluo_diffs + tr_probs_log;
        total_log_likelihoods = exp(total_log_likelihoods - logsumexp(total_log_likelihoods,1));
        
        %%% draw new samples
        option_array = cumsum(total_log_likelihoods);
        rand_array = repmat(rand(1,n_chains),nStates,1);
        sample_chains_slice(ind,:) = sum(rand_array > option_array) + 1;
    end
%     mcmcInfo.sample_chains(1,:,n) = randsample(state_options,n_chains,true,pi0);    
%     % simulate remaining time step
%     for t = 2:sql
%         % simulate transitions
%         A_array = cumsum(A(:,mcmcInfo.sample_chains(t-1,:,n)));
%         rand_array = repmat(rand(1,n_chains),nStates,1);
%         mcmcInfo.sample_chains(t,:,n) = sum(rand_array > A_array) + 1;
%     end    
end   
