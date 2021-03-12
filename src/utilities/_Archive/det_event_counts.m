function [n_init, n_steps, A_counts, posterior_prob_cell] = det_event_counts(...
    emissions_cell, fwd_probs, bkd_probs, A_curr, v_curr)

K = size(A_curr,1);
A_counts = zeros(K,K);
posterior_prob_cell = cell(size(emissions_cell));
% get counts of transitions
for e = 1:numel(emissions_cell) 
    % iterate through sequences
    init_vec = emissions_cell{e};
    fwd_array = fwd_probs{e};
    bkd_array = bkd_probs{e};  
    tr_p_mat = zeros(K^2,size(fwd_array,2)-1)-Inf;
    for k = 1:K
        for l = 1:K
            emission_rate = v_curr(l);
            akl = log(A_curr(l,k));                      
            % calculate 
            tr_p_mat((k-1)*K + l,:) = fwd_array(k,1:end-1) + akl + bkd_array(l,2:end)...
                + log(poisspdf(init_vec(2:end),emission_rate));            
        end        
    end
    tr_p_mat = tr_p_mat - logsumexp(tr_p_mat,1);
    % record
    A_counts = A_counts + reshape(sum(exp(tr_p_mat),2),K,K);
end

% get counts of emission events
n_init = zeros(K,1);
n_steps = zeros(K,1);
for e = 1:numel(emissions_cell)     
    init_vec = emissions_cell{e};
    fwd_array = fwd_probs{e};
    bkd_array = bkd_probs{e}; 
    prob_array = fwd_array + bkd_array;
    prob_array = exp(prob_array - logsumexp(prob_array,1));
    posterior_prob_cell{e} = prob_array;
    wt_array = init_vec .* prob_array;
    n_init = n_init + sum(wt_array,2);
    n_steps = n_steps + sum(prob_array,2);
end
