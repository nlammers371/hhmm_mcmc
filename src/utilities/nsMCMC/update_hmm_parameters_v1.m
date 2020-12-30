function mcmcInfo = update_hmm_parameters_v1(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    % update A
    A_counts = sum(mcmcInfo.transition_count_array,3)/mcmcInfo.n_chains;    
    mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts);
    mcmcInfo.A_inf_array(:,:,mcmcInfo.step) = mcmcInfo.A_curr;
    
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
            state_counts = convn(coeff_MS2,mcmcInfo.sample_chains==m,'full');            
            F_array(ind1:ind2,:,m) = state_counts(1:end-length(coeff_MS2)+1,n,:);                        
        end
    end
       
    F_array_long = reshape(F_array,[],3);
    y_vec = reshape(y_array,[],1);
    M = F_array_long'*F_array_long + realmin;    
    b = F_array_long'*y_vec;
    
    % calculate mean and variance
    v_mean = M\b;
    v_cov_mat = mcmcInfo.sigma^2 * inv(M);
    
    % sample
    mcmcInfo.v_curr = mvnrnd(v_mean, v_cov_mat)';
    mcmcInfo.v_inf_array(mcmcInfo.step,:) = mcmcInfo.v_curr;
    
    % get predicted fluorescence
    mcmcInfo = predict_fluo_full(mcmcInfo);
    
    % Update sigma
    a = numel(mcmcInfo.observed_fluo)/2;
    F_diff = reshape(mean(permute(mcmcInfo.sample_fluo,[1 3 2]),3) - mcmcInfo.observed_fluo,[],1);
    b = F_diff'*F_diff / 2;
%     y_vec = mcmcInfo.observed_fluo(:);    
%     F_diff = F_array_long * mcmcInfo.v_curr - y_vec;    
    mcmcInfo.sigma_curr = sqrt(1./gamrnd(a,1./b));%mcmcInfo.sigma;%sqrt(mean(F_diff.^2));%
    