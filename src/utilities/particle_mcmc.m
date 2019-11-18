function [v_cts_out, v_wts_out, A_count, pi_cell,prob_cell] = particle_mcmc_bkd(n_particles,...
    emissions_cell, A_curr, v_curr, pi0_curr)

K = size(A_curr,1);
A_tracker_array = NaN(K,K,numel(emissions_cell));
v_cts = NaN(K,1,numel(emissions_cell));
v_wts = NaN(K,1,numel(emissions_cell));
pi_cell = cell(1,numel(emissions_cell));
logL_cell = cell(1,numel(emissions_cell));
for e = 1:numel(emissions_cell)     
    init_vec = emissions_cell{e};
    T = numel(init_vec);
    % array to store particle states
    prob_array = zeros(n_particles,T);
    pt_array = zeros(n_particles,T);   
    occ_array = zeros(K,T);
    % index vector for sampling
    index_vec = 1:n_particles;
    %%%% initial state
    % sample initial state
    init_prop_vec = randsample(1:K,n_particles,true,pi0_curr);
    for k = 1:K
        occ_array(k,1) = sum(init_prop_vec==k);
    end
    % resample
    emit_probs = poisspdf(init_vec(1),v_curr(init_prop_vec));
    rs_indices = randsample(index_vec,numel(index_vec),true,emit_probs);
    pt_array(:,1) = init_prop_vec(rs_indices);
    % store probabilities
    prob_array(:,1) = emit_probs(rs_indices).*pi0_curr(pt_array(:,1))';
    A_tracker = zeros(K,K);
    %%%% iteration through time points
    for n = 2:T
        % simulate transitions
        A_array = cumsum(A_curr(:,pt_array(:,n-1)));
        rand_array = repmat(rand(1,n_particles),K,1);
        state_vec_prop = sum(rand_array > A_array) + 1;
        for k = 1:K
            occ_array(k,n) = sum(init_prop_vec==k);
        end
        % resample
        emit_probs = poisspdf(init_vec(n-1),v_curr(state_vec_prop));
        rs_indices = randsample(index_vec,numel(index_vec),true,emit_probs);        
        % track transitions
        rs_from = pt_array(rs_indices,n-1);
        rs_to = state_vec_prop(rs_indices)';   
        for i = 1:n_particles
            A_tracker(rs_to(i),rs_from(i)) = A_tracker(rs_to(i),rs_from(i)) + 1;
            % store probabilities
            prob_array(i,n) = A_curr(rs_to(i),rs_from(i))*emit_probs(rs_indices(i));
        end        
        pt_array(:,n) = rs_to;
    end     
    % record         
    pi_array = NaN(K,size(pt_array,2));
    for k = 1:K
        pi_array(k,:) = sum(pt_array == k);
    end
    pi_array = pi_array ./ size(pt_array,1);
    pi_cell{e} = pi_array;   
    prob_cell{e} = prob_array;
    %%% now check implied "A" and "v" values for consistency
    % v first   
    for k = 1:K
        v_cts(k,e) = sum((sum(pt_array==k) .* init_vec));
        v_wts(k,e) = sum(pt_array(:)==k);
    end
    % now A
    A_tracker_array(:,:,e) = A_tracker;        
end
% extract values
A_count = sum(A_tracker_array,3);
v_wts_out = sum(v_wts,2);
v_cts_out = sum(v_cts,2);