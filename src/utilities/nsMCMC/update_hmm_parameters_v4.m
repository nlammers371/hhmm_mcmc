function mcmcInfo = update_hmm_parameters_v4(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains_eff;
%     coeff_MS2 = mcmcInfo.coeff_MS2;        
    coeff_MS2_us = mcmcInfo.coeff_MS2_us;        
    us_factor = mcmcInfo.upsample_factor; 
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    mcmcInfo.update_flag = update_flag;
    if mcmcInfo.update_increment~=1
        update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    else
        update_index = mcmcInfo.step;
    end
    tr_counts = mcmcInfo.transition_count_array/us_factor;  
    pi0_counts = mcmcInfo.state_counts;
    if mcmcInfo.ensembleInferenceFlag
        tr_counts = repmat(mean(tr_counts,3),1,1,n_chains);
        pi0_counts = repmat(mean(pi0_counts,3),1,1,n_chains);
    end
    
%     ref_chain_ids = repelem(find(mcmcInfo.refChainVec),mcmcInfo.n_temps_per_chain);
    for n = 1:n_chains
        if n == 1 || ~mcmcInfo.ensembleInferenceFlag
            if ~mcmcInfo.rateSamplingFlag
                T = 1;%mcmcInfo.tempGradVec(n);
                A_chain = tr_counts(:,:,n).^(1/T); % NL: note that this technically should be applied to full distribution, including prior       
                A_samp = sample_A_dirichlet(mcmcInfo.A_alpha(:,:,n), A_chain);    
                mcmcInfo.A_curr(:,:,n) = A_samp;
            else
                if nStates == 2
                    kon_tr = tr_counts(2,1,n)+mcmcInfo.alpha_kon;
                    kon_total = pi0_counts(1,1,n)+mcmcInfo.beta_kon;
                    koff_tr = tr_counts(1,2,n)+mcmcInfo.alpha_koff;
                    koff_total = pi0_counts(1,2,n)+mcmcInfo.beta_koff;
                    % sample from gamma 
                    kon = gamrnd(kon_tr,1/kon_total);
                    koff = gamrnd(koff_tr,1/koff_total);
                    % generate new rate matrix
                    Q_init = [-kon koff;  kon -koff];
                    mcmcInfo.Q_curr(:,:,n) = Q_init/mcmcInfo.tres*mcmcInfo.upsample_factor;
                    mcmcInfo.A_curr(:,:,n) = expm(Q_init);%eye(2) + Q_init + Q_init^2/2 + Q_init^3/6 + Q_init^4/24;
                    
                elseif nStates == 3
                  
                    % calculate for high and low states first
                    kon_tr = tr_counts(2,1,n)+2*mcmcInfo.alpha_kon;
                    kon_total = pi0_counts(1,1,n)+mcmcInfo.beta_kon;
                    koff_tr = tr_counts(2,3,n)+2*mcmcInfo.alpha_koff;
                    koff_total = pi0_counts(1,3,n)+mcmcInfo.beta_koff;
                    
                    % sample from gamma 
                    kon = gamrnd(kon_tr,1/kon_total);
                    koff = gamrnd(koff_tr,1/koff_total);
                    
                    % now middle rates
                    kon_tr2 = tr_counts(3,2,n)+mcmcInfo.alpha_kon;
                    koff_tr2 = tr_counts(1,2,n)+mcmcInfo.alpha_koff;
                    mid_total = pi0_counts(1,2,n)+mcmcInfo.beta_koff+mcmcInfo.beta_kon;
          
                    % sample from gamma. Note that we sample combined rate
                    % out
                    k_mid = gamrnd(kon_tr2+koff_tr2,1/mid_total);
                    k_ratio = kon_tr2 / (kon_tr2 + koff_tr2);
                    kon2 = k_mid*k_ratio;
                    koff2 = k_mid*(1-k_ratio);
                    
                    % make rate matrix
                    Q_init = [-kon koff2 0; kon -k_mid koff; 0 kon2 -koff];
                    mcmcInfo.Q_curr(:,:,n) = Q_init/mcmcInfo.tres*mcmcInfo.upsample_factor;
                    if mcmcInfo.adjustSamplingFlag
                        mcmcInfo.A_curr(:,:,n) = expm(Q_init);%eye(3) + Q_init + Q_init^2/2 + Q_init^3/6 + Q_init^4/24;                    
                    else
                        mcmcInfo.A_curr(:,:,n) = eye(3) + Q_init;
                    end
                    mcmcInfo.A_curr(:,:,n) = mcmcInfo.A_curr(:,:,n) ./ sum(mcmcInfo.A_curr(:,:,n),1);
                else
                    error('Rate sampling not supported for nStates>3');
                end
            end
 
            % update pi0      
            pi0_ct = pi0_counts(1,:,n) + sum(mcmcInfo.A_alpha(:,:,n),1);
            mcmcInfo.pi0_curr(n,:) = drchrnd(pi0_ct,1);
        else
            mcmcInfo.A_curr(:,:,n) = mcmcInfo.A_curr(:,:,1);
            mcmcInfo.pi0_curr(n,:) = mcmcInfo.pi0_curr(1,:);
        end

        % check that pi0 values are pos
        if update_flag
            mcmcInfo.A_inf_array(:,:,update_index,n) = mcmcInfo.A_curr(:,:,n);            
            mcmcInfo.pi0_inf_array(update_index,:,n) = mcmcInfo.pi0_curr(n,:);
            if mcmcInfo.rateSamplingFlag
                mcmcInfo.Q_inf_array(:,:,update_index,n) = mcmcInfo.Q_curr(:,:,n);
            end
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
                state_counts = convn(coeff_MS2_us(:,c),mcmcInfo.sample_chains(:,c,n)==m,'full')/us_factor;            
                state_counts = state_counts(1:end-size(coeff_MS2_us,1)+1,:);
                state_counts = state_counts(us_factor:us_factor:end,:);
                F_array(ind1:ind2,m,c) = state_counts;                        
            end
        end
    end  
    if mcmcInfo.ensembleInferenceFlag   
        F_array = repmat(mean(F_array,3),1,1,n_chains);
    end
    
    swap_indices = NaN(n_chains,nStates);
    
    for c = 1:n_chains   
        if c == 1 || ~mcmcInfo.ensembleInferenceFlag
            T = 1;%mcmcInfo.tempGradVec(c);
            M = ((F_array(:,:,c)'*F_array(:,:,c))) + 1e-1;    
            b = ((F_array(:,:,c)'*y_array(:,c)));

            % calculate mean and variance
            v_lsq = M\b;
            v_mean = (M + mcmcInfo.M0)^-1 * (mcmcInfo.M0*mcmcInfo.v0(c,:)' + M*v_lsq);
            v_cov_mat = T * inv(mcmcInfo.sigma_curr(c)^-2 * M +  mcmcInfo.sigma_curr(c)^-2 * inv(mcmcInfo.M0));

            % sample
            try
            [mcmcInfo.v_curr(c,:), swap_indices(c,:)] = sort(mvnrnd(v_mean, v_cov_mat)'); 
            catch
                error('wtf')
%                 v_cov_mat = T * inv(mcmcInfo.sigma_curr(c)^-2 *inv(mcmcInfo.M0));
%                 mcmcInfo.v_curr(c,:) = mvnrnd(v_mean, v_cov_mat)'; 
%                 mcmcInfo.err_flag_vec(c) = 1;
            end
        else
            mcmcInfo.v_curr(c,:) = mcmcInfo.v_curr(1,:);
        end
        if update_flag
            mcmcInfo.v_inf_array(update_index,:,c) = mcmcInfo.v_curr(c,:);   
        end
    end
    
    % deal with label switches
    for c = 1:n_chains
        si = swap_indices(c,:);
        mcmcInfo.A_curr(:,:,c) = mcmcInfo.A_curr(si,si,c);        
        mcmcInfo.pi0_curr(c,:) = mcmcInfo.pi0_curr(c,si);
        
        sc = mcmcInfo.sample_chains(:,c,:);
        sc_new = sc;
        for k = 1:nStates
            sc_new(sc==k) = si(k);
        end
        mcmcInfo.sample_chains(:,c,:) = sc_new; 
        
        if update_flag
            mcmcInfo.pi0_inf_array(update_index,:,c) = mcmcInfo.pi0_inf_array(update_index,si,c);
            mcmcInfo.v_inf_array(update_index,:,c) = mcmcInfo.v_inf_array(update_index,si,c);
            mcmcInfo.A_inf_array(:,:,update_index,c) = mcmcInfo.A_inf_array(si,si,update_index,c);
        end
    end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%
    
    % get predicted fluorescence (using new v values)
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);
    
    if mcmcInfo.ensembleInferenceFlag
        mcmcInfo.sample_fluo = repmat(mean(mcmcInfo.sample_fluo,2),1,n_chains);
    end
    
    % Update sigma
    for c = 1:n_chains
        if c == 1 || ~mcmcInfo.ensembleInferenceFlag
            T = 1;%mcmcInfo.tempGradVec(c);
            % see: https://discdown.org/flexregression/bayesreg.html                        
            if ~mcmcInfo.bootstrapFlag
                a = (numel(mcmcInfo.observed_fluo)/2 + mcmcInfo.a0)./T;    
                F_diff = reshape(permute(mcmcInfo.sample_fluo(us_factor:us_factor:end,c,:),[1 3 2]) - mcmcInfo.observed_fluo,[],1);       
            else
                F_diff = reshape(mcmcInfo.sample_fluo(us_factor:us_factor:end,c,:) - mcmcInfo.observed_fluo(:,c,:),[],1);       
                a = (numel(mcmcInfo.observed_fluo(:,1,:))/2 + mcmcInfo.a0)./T;    
            end
            b_prior_piece = mcmcInfo.b0 + (mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))*inv(mcmcInfo.M0)*(mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))';
            b = (F_diff'*F_diff / 2 + b_prior_piece);

            mcmcInfo.sigma_curr(c) = sqrt(1./gamrnd(a,1./b));%mcmcInfo.trueParams.sigma
        else
            mcmcInfo.sigma_curr(c) = mcmcInfo.sigma_curr(1);
        end
        if update_flag
            mcmcInfo.sigma_inf_array(update_index,c) = mcmcInfo.sigma_curr(c);
        end
    end   
    
    if update_flag && mcmcInfo.save_trace_results
        mcmcInfo.sample_fluo_inf_array(:,:,:,update_index) = mcmcInfo.sample_fluo;
        mcmcInfo.sample_states_inf_array(:,:,:,update_index) = mcmcInfo.sample_chains;
    end