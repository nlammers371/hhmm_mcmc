function [n_init, n_steps, A_counts, posterior_prob_cell] = particle_event_counts(...
    emissions_cell, fwd_probs, bkd_probs, A_curr, v_curr)

K = size(A_curr,1);
A_counts = zeros(K,K);
posterior_prob_cell = cell(size(emissions_cell));
% get counts of transitions
for e = 1:numel(emissions_cell) 
    % iterate through sequences
    init_vec = emissions_cell{e};
    fwd_array = log(fwd_probs{e});
    bkd_array = log(bkd_probs{e});  
    tr_p_mat = zeros(K,K,size(fwd_array,2)-1);
    for k = 1:K
        for l = 2:K
            emission_rate = v_curr(l);
            akl = log(A_curr(l,k));                      
            % calculate 
            tr_p_mat(l,k,:) = fwd_array(k,1:end-1) + akl + bkd_array(l,2:end)...
                + log(poisspdf(init_vec(2:end),emission_rate));            
        end        
    end
    tr_p_mat = tr_p_mat ./ sum(sum(tr_p_mat,1),2);
    % record
    A_counts = A_counts + sum(tr_p_mat,3);
end

% get counts of emission events
n_init = zeros(K,1);
n_steps = zeros(K,1);
for e = 1:numel(emissions_cell)     
    init_vec = emissions_cell{e};
    fwd_array = fwd_probs{e}+1e-6;
    bkd_array = bkd_probs{e}+1e-6; 
    prob_array = fwd_array.*bkd_array;
    prob_array = prob_array ./ sum(prob_array);
    posterior_prob_cell{e} = prob_array;
    wt_array = init_vec .* prob_array;
    n_init = n_init + sum(wt_array,2);
    n_steps = n_steps + sum(prob_array,2);
end
