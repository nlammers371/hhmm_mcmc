function mcmcInfo = update_hmm_parameters_v3(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    seq_length_obs = mcmcInfo.seq_length_obs;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains_eff;
    coeff_MS2 = mcmcInfo.coeff_MS2;        
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    mcmcInfo.update_flag = update_flag;
    if mcmcInfo.update_increment~=1
        update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    else
        update_index = mcmcInfo.step;
    end    
    A_counts = mcmcInfo.transition_count_array;  
%     ref_chain_ids = repelem(find(mcmcInfo.refChainVec),mcmcInfo.n_temps_per_chain);
    for n = 1:n_chains
        T = 1;%mcmcInfo.tempGradVec(n);
        A_chain = A_counts(:,:,n).^(1/T); % NL: note that this technically should be applied to full distribution, including prior       
        A_samp = sample_A_dirichlet(mcmcInfo.A_alpha(:,:,n), A_chain);    
        mcmcInfo.A_curr(:,:,n) = A_samp;
        
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
    y_array = NaN(seq_length_obs*n_traces,n_chains);    
        
    for c = 1:n_chains
        for n = 1:n_traces        
            ind1 = (n-1)*seq_length+1;
            ind2 = n*seq_length;
            ind1_obs = (n-1)*seq_length_obs+1;
            ind2_obs = n*seq_length_obs;
            % record observed fluo
            y_array(ind1_obs:ind2_obs,c) = mcmcInfo.observed_fluo(:,n);
            for m = 1:nStates
                % record counts
                state_counts = convn(coeff_MS2(:,c),mcmcInfo.sample_chains(:,c,n)==m,'full');            
                F_array(ind1:ind2,m,c) = state_counts(1:end-size(coeff_MS2,1)+1,:);                        
            end
        end
    end  
%     if mcmcInfo.upsample_factor > 1
%         conv_kernel = ones(mcmcInfo.upsample_factor,1);
%         F_array_conv = convn(conv_kernel,F_array,'full');
%         F_array = F_array_conv(mcmcInfo.upsample_factor:mcmcInfo.upsample_factor:end,:,:);
%     end
    F_array = F_array(mcmcInfo.upsample_factor:mcmcInfo.upsample_factor:end,:,:);
    
    for c = 1:n_chains        
        T = 1;%mcmcInfo.tempGradVec(c);
        M = ((F_array(:,:,c)'*F_array(:,:,c))) + 1e-1;    
        b = ((F_array(:,:,c)'*y_array(:,c)));
                 
        % calculate mean and variance
        v_lsq = M\b;
        v_mean = (M + mcmcInfo.M0)^-1 * (mcmcInfo.M0*mcmcInfo.v0(c,:)' + M*v_lsq);
        v_cov_mat = T * inv(mcmcInfo.sigma_curr(c)^-2 * M +  mcmcInfo.sigma_curr(c)^-2 *inv(mcmcInfo.M0));
        
        % sample
        mcmcInfo.v_curr(c,:) = mvnrnd(v_mean, v_cov_mat)'; 

        if update_flag
            mcmcInfo.v_inf_array(update_index,:,c) = mcmcInfo.v_curr(c,:);   
        end
    end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%

    % get predicted fluorescence (using new v values)
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);
    
    % Update sigma
    for c = 1:n_chains
        T = 1;%mcmcInfo.tempGradVec(c);
        % see: https://discdown.org/flexregression/bayesreg.html
        a = (numel(mcmcInfo.observed_fluo)/2 + mcmcInfo.a0)./T;    
        
        F_diff = reshape(permute(mcmcInfo.sample_fluo(:,c,:),[1 3 2]) - mcmcInfo.observed_fluo,[],1);       
        b_prior_piece = mcmcInfo.b0 + (mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))*inv(mcmcInfo.M0)*(mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))';
        b = (F_diff'*F_diff / 2 + b_prior_piece);
        
        mcmcInfo.sigma_curr(c) = sqrt(1./gamrnd(a,1./b));%mcmcInfo.trueParams.sigma
        if update_flag
            mcmcInfo.sigma_inf_array(update_index,c) = mcmcInfo.sigma_curr(c);
        end
    end   