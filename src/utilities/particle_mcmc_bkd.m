function [pi_cell] = particle_mcmc_bkd(n_particles,emissions_cell, A_curr, v_curr, pi0_curr)

% take transpose of transition prob array
A_curr_bkd = A_curr' ./sum(A_curr');
% initialize shit
K = size(A_curr_bkd,1);
pi_cell = cell(1,numel(emissions_cell));
% iterate through traces
for e = 1:numel(emissions_cell)     
    init_vec = emissions_cell{e};
    T = numel(init_vec);
    % array to store particle states    
    pt_array = zeros(n_particles,T);   
    occ_array = zeros(K,T);
    % index vector for sampling
    index_vec = 1:n_particles;
    %%%% initial state
    % sample initial state
    init_prop_vec = randsample(1:K,n_particles,true,pi0_curr);        
    % record state occupancy
    pt_array(:,T) = init_prop_vec;
    for k = 1:K
        occ_array(k,1) = sum(pt_array(:,T)==k);
    end    
    %%%% iteration through time points
    for n = fliplr(1:T-1)
        % simulate transitions
        A_array = cumsum(A_curr_bkd(:,pt_array(:,n+1)));
        rand_array = repmat(rand(1,n_particles),K,1);
        state_vec_prop = sum(rand_array > A_array) + 1;        
        % resample
        try
            emit_probs = poisspdf(init_vec(n+1),v_curr(state_vec_prop));
        catch
            error('sigh')
        end
        rs_indices = randsample(index_vec,numel(index_vec),true,emit_probs);         
        pt_array(:,n) = state_vec_prop(rs_indices);
        for k = 1:K
            occ_array(k,n) = sum(pt_array(:,n)==k);
        end
    end     
    % record         
    pi_array = NaN(K,size(pt_array,2));
    for k = 1:K
        pi_array(k,:) = sum(pt_array == k);
    end
    pi_array = pi_array ./ size(pt_array,1);
    pi_cell{e} = pi_array;       
end