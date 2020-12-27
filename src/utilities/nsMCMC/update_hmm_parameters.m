function mcmcInfo = update_hmm_parameters(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
%     n_chains = mcmcInfo.n_chains;
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    % update A
    A_counts = sum(mcmcInfo.transition_count_array,3)/mcmcInfo.n_chains;    
    mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts);
    mcmcInfo.A_inf_array(:,:,mcmcInfo.step) = mcmcInfo.A_curr;
    
    %%% use ML formula to sample emissions %%%
    
    % generate F count arrays
    F_array = zeros(seq_length,nStates,n_traces);    
        
    for m = 1:nStates
        state_counts = convn(coeff_MS2,mcmcInfo.sample_chains==m,'full');
        F_array(:,m,:) = mean(state_counts(1:end-length(coeff_MS2)+1,:,:),2);
    end
    
    % generate M and b arrays    
    M = zeros(nStates);
    b = zeros(nStates,1);
    y_array = mcmcInfo.observed_fluo;
    for m = 1:nStates
        for n = 1:nStates
            M(m,n) = sum(sum(F_array(:,m,:).*F_array(:,n,:)),3);
            b(m) = sum(sum(permute(F_array(:,m,:),[1 3 2]).*y_array),2);
        end
    end
    % update
    mcmcInfo.v_curr = M\b;
    mcmcInfo.v_inf_array(mcmcInfo.step,:) = mcmcInfo.v_curr;