function mcmcInfo = update_hmm_parameters_PF(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains_eff;
%     coeff_MS2 = mcmcInfo.coeff_MS2;        
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    mcmcInfo.update_flag = update_flag;
    if mcmcInfo.update_increment~=1
        update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    else
        update_index = mcmcInfo.step;
    end
    A_counts = mcmcInfo.transition_count_array;  
    pi0_counts = mcmcInfo.pi0_counts;
                
    T = 1; 
    A_chain = A_counts.^(1/T); % NL: note that this technically should be applied to full distribution, including prior       
    A_samp = sample_A_dirichlet(mcmcInfo.A_alpha(:,:,1), A_chain);    
    mcmcInfo.A_curr(:,:,:) = repmat(A_samp,1,1,n_chains);

    % update pi0      
    pi0_ct = pi0_counts' + sum(mcmcInfo.A_alpha(:,:,1),1);
    mcmcInfo.pi0_curr(:,:) = repmat(drchrnd(pi0_ct,1),n_chains,1);
       
    if update_flag
        mcmcInfo.A_inf_array(:,:,update_index,:) = mcmcInfo.A_curr;
        mcmcInfo.pi0_inf_array(update_index,:,:) = permute(mcmcInfo.pi0_curr,[3 2 1]);
    end                    
    
    %% %%%%%%%%%%%%% update emission vector(V) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V   
    
    % generate F count arrays
    F_array_deep = mcmcInfo.F_array;
    F_array = zeros(seq_length*n_traces,nStates);        
    y_array = NaN(seq_length*n_traces,1);    
            
    for n = 1:n_traces        
        ind1 = (n-1)*seq_length+1;
        ind2 = n*seq_length;
        
        % record observed fluo                    
        y_array(ind1:ind2,1) = mcmcInfo.observed_fluo(:,n);                    
        F_array(ind1:ind2,:) = F_array_deep(:,:,n);                        
        
    end
    
    if mcmcInfo.ensembleInferenceFlag
        if mcmcInfo.bootstrapFlag
            error('incompatible inference options')
        end
        F_array = repmat(mean(F_array,3),1,1,n_chains);
    end

    T = 1;%mcmcInfo.tempGradVec(c);
    M = ((F_array'*F_array)) + 1e-1;    
    b = ((F_array'*y_array));

    % calculate mean and variance
    v_lsq = M\b;
    v_mean = (M + mcmcInfo.M0)^-1 * (mcmcInfo.M0*mcmcInfo.v0(1,:)' + M*v_lsq);
    v_cov_mat = T * inv(mcmcInfo.sigma_curr(1)^-2 * M +  mcmcInfo.sigma_curr(1)^-2 * inv(mcmcInfo.M0));

    % sample
    try
        mcmcInfo.v_curr(:,:) = repmat(mvnrnd(v_mean, v_cov_mat),n_chains,1); 
    catch
        error('wtf')
    end
    
    if update_flag
        mcmcInfo.v_inf_array(update_index,:,:) = permute(mcmcInfo.v_curr(:,:),[3 2 1]);   
    end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%

    % get predicted fluorescence (using new v values)
%     mcmcInfo = predict_fluo_full_v3(mcmcInfo);
    
%     if mcmcInfo.ensembleInferenceFlag
%         mcmcInfo.sample_fluo = repmat(mean(mcmcInfo.sample_fluo,2),1,n_chains);
%     end
    
    % Update sigma    
    T = 1;%mcmcInfo.tempGradVec(c);
    % see: https://discdown.org/flexregression/bayesreg.html                            
    a = (numel(mcmcInfo.observed_fluo)/2 + mcmcInfo.a0)./T;    
    F_diff = (mcmcInfo.sample_fluo - permute(mcmcInfo.observed_fluo,[1,3,2])).^2;       
    F_diff_avg = sum(exp(mcmcInfo.sample_logL_bkd).*F_diff,2)./sum(exp(mcmcInfo.sample_logL_bkd),2);
    b_prior_piece = mcmcInfo.b0 + (mcmcInfo.v_curr(1,:)-mcmcInfo.v0(1,:))*inv(mcmcInfo.M0)*(mcmcInfo.v_curr(1,:)-mcmcInfo.v0(1,:))';
    b = (sum(F_diff_avg(:)) / 2 + b_prior_piece);

    mcmcInfo.sigma_curr = repmat(sqrt(1./gamrnd(a,1./b)),n_chains,1);%mcmcInfo.trueParams.sigma
    if update_flag
        mcmcInfo.sigma_inf_array(update_index,:) = mcmcInfo.sigma_curr(1);
    end
    
    if update_flag && mcmcInfo.save_trace_results
        mcmcInfo.sample_fluo_inf_array(:,:,:,update_index) = mcmcInfo.sample_fluo;
        mcmcInfo.sample_states_inf_array(:,:,:,update_index) = mcmcInfo.sample_chains;
    end