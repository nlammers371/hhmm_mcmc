function [pi_cell,logL_seq,logL_tot] = particle_mcmc_fwd(n_particles,...
    emissions_cell, A_curr, v_curr, pi0_curr)

K = size(A_curr,1);
A_tracker_array = NaN(K,K,numel(emissions_cell));
v_cts = NaN(K,1,numel(emissions_cell));
v_wts = NaN(K,1,numel(emissions_cell));
pi_cell = cell(1,numel(emissions_cell));
logL_seq = NaN(1,numel(emissions_cell));
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
    % resample
    emit_probs = poisspdf(init_vec(1),v_curr(init_prop_vec));
    rs_indices = randsample(index_vec,numel(index_vec),true,emit_probs);
    pt_array(:,1) = init_prop_vec(rs_indices);
    for k = 1:K
        occ_array(k,1) = sum(pt_array(:,1)==k);
    end
    % store probabilities
    prob_array(:,1) = emit_probs(rs_indices).*pi0_curr(pt_array(:,1))';
    A_tracker = zeros(K,K);
    %%%% iteration through time points
    for n = 2:T
        % simulate transitions
        A_array = cumsum(A_curr(:,pt_array(:,n-1)));
        rand_array = repmat(rand(1,n_particles),K,1);
        state_vec_prop = sum(rand_array > A_array) + 1;        
        % resample
        emit_probs = poisspdf(init_vec(n),v_curr(state_vec_prop));
        rs_indices = randsample(index_vec,numel(index_vec),true,emit_probs);        
        % track transitions
        rs_from = pt_array(rs_indices,n-1);
        rs_to = state_vec_prop(rs_indices)'; 
        lin_ind = sub2ind(size(A_curr), rs_to, rs_from);     
        % record
        prob_array(:,n) = A_curr(lin_ind).*emit_probs(rs_indices)';            
        pt_array(:,n) = rs_to;
        % get occupancy
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
    logL_seq(e) = sum(mean(log(prob_array)));           
end
logL_tot = sum(logL_seq);