function [bkd_probs] = bkd_algorithm(seq_cell, A_curr, v_curr, pi0_curr)
    % get N states
    K = size(A_curr,1);
    % get log versions of matrices
    A_log = log(A_curr);    
    pi0_log = log(pi0_curr);
    % get N seqs 
    n_seq = numel(seq_cell);
    % initialize arrays
    bkd_probs = cell(size(seq_cell));    
    for s = 1:n_seq        
        init_seq = seq_cell{s};
        T = numel(init_seq);
        % precalculate array of emission probabilities
        init_probs_log = log(poisspdf(repmat(init_seq,K,1),repmat(v_curr',1,T)));
        % initialize bkd array                
        bkd_array = NaN(K,T);
        bkd_array(:,T) = pi0_log;
        % iterate        
        for t = fliplr(1:T-1)                        
            post = repmat(bkd_array(:,t+1),1,K);
            bkd_array(:,t) = logsumexp(post + A_log + repmat(init_probs_log(:,t+1),1,K),1);            
        end
        bkd_probs{s} = bkd_array;
    end    
end
