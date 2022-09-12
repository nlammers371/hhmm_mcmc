function mcmcInfo = gibbs_update_A(mcmcInfo)

    us_factor = mcmcInfo.upsample_factor;
    n_chains = mcmcInfo.n_chains;
    nStates = mcmcInfo.nStates;
    update_flag = mcmcInfo.update_flag;
    
    if mcmcInfo.update_increment~=1
        update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    else
        update_index = mcmcInfo.step;
    end
    
    tr_counts = mcmcInfo.transition_count_array/us_factor;  
    pi0_counts = mcmcInfo.state_counts/us_factor;
    for n = 1:n_chains

        if ~mcmcInfo.rateSamplingFlag                
            A_chain = tr_counts(:,:,n); 
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
                mcmcInfo.A_curr(:,:,n) = eye(2) + Q_init;%expm(Q_init);%eye(2) + Q_init + Q_init^2/2 + Q_init^3/6 + Q_init^4/24;

            elseif nStates == 3

                % calculate for high and low states first
                kon_tr = tr_counts(2,1,n)+2*mcmcInfo.alpha_kon;
                kon_total = pi0_counts(1,1,n)+2*mcmcInfo.beta_kon;
                koff_tr = tr_counts(2,3,n)+2*mcmcInfo.alpha_koff;
                koff_total = pi0_counts(1,3,n)+2*mcmcInfo.beta_koff;

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
                if ~mcmcInfo.adjustSamplingFlag
                    mcmcInfo.A_curr(:,:,n) = generate_A_matrix(Q_init);                
                else
                    mcmcInfo.A_curr(:,:,n) = expm(Q_init);                
                end

            else
                error('Rate sampling not supported for nStates>3');
            end
        end

        % update pi0      
        pi0_ct = pi0_counts(1,:,n) + sum(mcmcInfo.A_alpha(:,:,n),1);
        mcmcInfo.pi0_curr(n,:) = drchrnd(pi0_ct,1);

        % check that pi0 values are pos
        if update_flag
            mcmcInfo.A_inf_array(:,:,update_index,n) = mcmcInfo.A_curr(:,:,n);            
            mcmcInfo.pi0_inf_array(update_index,:,n) = mcmcInfo.pi0_curr(n,:);
            if mcmcInfo.rateSamplingFlag
                mcmcInfo.Q_inf_array(:,:,update_index,n) = mcmcInfo.Q_curr(:,:,n);
            end
        end
    end