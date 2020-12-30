function mcmcInfo = update_hmm_parameters_gibbs(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    coeff_MS2 = mcmcInfo.coeff_MS2;
    

    % update A
    A_counts = mean(mcmcInfo.transition_count_array,3)/mcmcInfo.n_chains;    
    A_ml = mcmcInfo.A_alpha + A_counts;
    A_ml = A_ml ./ sum(A_ml);
    mcmcInfo.A_inf_array(:,:,mcmcInfo.step) = A_ml;      
    if mcmcInfo.par_chain_flag
        mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts, n_chains);        
    else
        mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts, 1);          
    end
    
    % update V
    
    % generate F count arrays
    F_array = zeros(seq_length*n_chains,n_traces,nStates);    
    y_array = NaN(seq_length*n_chains,n_traces);    
    
    for n = 1:n_chains
        ind1 = (n-1)*seq_length+1;
        ind2 = n*seq_length;
        % record observed fluo
        y_array(ind1:ind2,:) = mcmcInfo.observed_fluo;
        for m = 1:nStates
            % record counts
            state_counts = convn(coeff_MS2,mcmcInfo.sample_chains(:,n,:)==m,'full');            
            F_array(ind1:ind2,:,m) = state_counts(1:end-length(coeff_MS2)+1,1,:);                        
        end
    end
       
    F_array_long = reshape(F_array,[],3);
    y_vec = reshape(y_array,[],1);
    M = F_array_long'*F_array_long + realmin;    
    b = F_array_long'*y_vec;
    
    % calculate mean and variance
    v_array_prev = mcmcInfo.v_curr;
    v_ml = M\b;
    v_cov_mat = mcmcInfo.sigma_inf_array(mcmcInfo.step-1)^2 * inv(M); % NL: I'm using previous ML sigma value for this    
    
    % sample
    mcmcInfo.v_inf_array(mcmcInfo.step,:) = v_ml;
    if mcmcInfo.par_chain_flag
        mcmcInfo.v_curr = mvnrnd(v_ml, v_cov_mat,mcmcInfo.n_chains);
    else
        mcmcInfo.v_curr = mvnrnd(v_ml, v_cov_mat,1);
    end
        
    % get predicted fluorescence
    predict_fluo_full(mcmcInfo, v_array_prev);
    
    % Update sigma
    a = numel(mcmcInfo.observed_fluo)/2;
    F_diff = reshape(mean(permute(mcmcInfo.sample_fluo,[1 3 2]),3) - mcmcInfo.observed_fluo,[],1);
    b = F_diff'*F_diff / 2;
  
    if mcmcInfo.par_chain_flag
        mcmcInfo.sigma_curr = sqrt(1./gamrnd(a,1./b,1,mcmcInfo.n_chains));
    else
        mcmcInfo.sigma_curr = sqrt(1./gamrnd(a,1./b));%mcmcInfo.sigma;%sqrt(mean(F_diff.^2));%
    end
    mcmcInfo.sigma_inf_array(mcmcInfo.step) = sqrt(mean(F_diff.^2));