function [fwd_probs, logL_seq, logL_tot] = fwd_algorithm(seq_cell, A_curr, v_curr, pi0_curr)
    % get N states
    K = size(A_curr,1);
    % get log versions of matrices
    A_log = log(A_curr);    
    pi0_log = log(pi0_curr);
    % get N seqs 
    n_seq = numel(seq_cell);
    % initialize arrays
    fwd_probs = cell(size(seq_cell));
    logL_seq = NaN(size(seq_cell));
    for s = 1:n_seq        
        init_seq = seq_cell{s};
        T = numel(init_seq);
        % precalculate array of emission probabilities
        init_probs_log = log(poisspdf(repmat(init_seq,K,1),repmat(v_curr',1,T)));
        % inistialize fwd array
        fwd_array = NaN(K,T);        
        % iterate
        prev = repmat(pi0_log',K,1);                
        for t = 1:T            
            fwd_array(:,t) = logsumexp(prev + A_log,2) + init_probs_log(:,t);
            prev = repmat(fwd_array(:,t)',K,1);
        end       
        fwd_probs{s} = fwd_array;
        logL_seq(s) = logsumexp(fwd_array(:,end),1);        
    end    
    logL_tot = sum(logL_seq);
end
