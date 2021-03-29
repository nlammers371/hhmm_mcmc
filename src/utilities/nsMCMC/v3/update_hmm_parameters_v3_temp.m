function mcmcInfo = update_hmm_parameters_v3_temp(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    coeff_MS2 = mcmcInfo.coeff_MS2;        
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    A_counts = mcmcInfo.transition_count_array;  
    for n = 1:n_chains
        mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet(mcmcInfo.A_alpha(:,:,n), A_counts(:,:,n));    
        
        % update pi0    
        [V, D] = eig(mcmcInfo.A_curr(:,:,n));
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));    
        
        if update_flag
            mcmcInfo.A_inf_array(:,:,update_index,n) = mcmcInfo.A_curr(:,:,n);
            mcmcInfo.pi0_inf_array(update_index,:,n) = mcmcInfo.pi0_curr(n,:);
        end
    end
                
    
    %% %%%%%%%%%%%%% update emission vector(V) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V   
    
    % generate F count arrays
    F_array = zeros(seq_length*n_traces,nStates,n_chains);    
    y_array = NaN(seq_length*n_traces,n_chains);    
    
    for c = 1:n_chains
        for n = 1:n_traces        
            ind1 = (n-1)*seq_length+1;
            ind2 = n*seq_length;
            % record observed fluo
            y_array(ind1:ind2,c) = mcmcInfo.observed_fluo(:,n);
            for m = 1:nStates
                % record counts
                state_counts = convn(coeff_MS2,mcmcInfo.sample_chains(:,c,n)==m,'full');            
                F_array(ind1:ind2,m,c) = reshape(state_counts(1:end-length(coeff_MS2)+1,:),[],1);                        
            end
        end
    end  
    
    for c = 1:n_chains
        % update for each chain       
        M = (F_array(:,:,c)'*F_array(:,:,c) + 1e-4);    
        b = (F_array(:,:,c)'*y_array(:,c));

        % calculate mean and variance
        v_mean = M\b;
        v_cov_mat = mcmcInfo.sigma_curr(c)^2 * inv(M);

        % sample
        mcmcInfo.v_curr(c,:) = mvnrnd(v_mean, v_cov_mat)';             

        if update_flag
            mcmcInfo.v_inf_array(update_index,:,c) = mcmcInfo.v_curr(c,:);   
        end
    end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);
       
    % Update sigma
    for c = 1:n_chains
        a = numel(mcmcInfo.observed_fluo)/2;
%         mean_fluo_est = permute(mean(mcmcInfo.sample_fluo,2),[1 3 2]);
%         F_diff = reshape(mean_fluo_est - mcmcInfo.observed_fluo,[],1);
        F_diff = reshape(permute(mcmcInfo.sample_fluo(:,c,:),[1 3 2]) - mcmcInfo.observed_fluo,[],1);
        b = F_diff'*F_diff / 2;
%         if mcmcInfo.step>100
%           error('wtf')
%         end
        % draw sample        
        if mcmcInfo.step < 350
            sigma_est = sqrt(1./gamrnd(a,1./b));
            mcmcInfo.sigma_curr(c) = sigma_est*(1+1*100./(100+mcmcInfo.step));%mcmcInfo.trueParams.sigma
        else
            mcmcInfo.sigma_curr(c) = sqrt(1./gamrnd(a,1./b));%mcmcInfo.trueParams.sigma
        end

        if update_flag
            mcmcInfo.sigma_inf_array(update_index,c) = mcmcInfo.sigma_curr(c);
        end
    end   