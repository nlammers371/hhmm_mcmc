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
        % initialize bkd array
        bkd_array = NaN(K,T);
        % iterate
        post = repmat(pi0_log',K,1);
        for t = fliplr(1:T)            
            init_probs_log = log(poisspdf(init_seq(t),v_curr))';
            bkd_array(:,t) = logsumexp(post + A_log,2) + init_probs_log;
            post = repmat(bkd_array(:,t)',K,1);
        end
        bkd_probs{s} = bkd_array;
    end    
end
