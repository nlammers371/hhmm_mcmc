function mcmcInfo = adjust_tr_counts_v2(mcmcInfo)

    nStates = mcmcInfo.nStates;    
    n_chains = mcmcInfo.n_chains;    
    us_factor = mcmcInfo.upsample_factor;
    
    % Calculate expected state visits and state transitions
    us_factor_fine = 100;
    us_tres_fine = mcmcInfo.tres/us_factor_fine/us_factor;
    Q = mcmcInfo.Q_curr;
    
    iter_vec_fw = reshape(1:us_factor_fine,1,1,[]);
    iter_vec_bk = reshape(us_factor_fine:-1:1,1,1,[]);

    % initialize new array to store state probs
    sample_chains_prob = cell(1,nStates);
    for s = 1:nStates
        sample_chains_prob{s} = zeros(size(mcmcInfo.sample_chains));
    end
    
    for n = 1:n_chains
        
      % extract transition prob matrix
        A_temp = generate_A_matrix(Q(:,:,n)*us_tres_fine);
        
        % initialize helper arrays for calculation\
        A_array_fw = NaN(nStates,nStates,us_factor_fine-1);
        A_array_bk = NaN(nStates,nStates,us_factor_fine-1);
        
        % generate to and from reference arrays
        from_array = mcmcInfo.sample_chains(1:end-1,n,:);
        to_array = mcmcInfo.sample_chains(2:end,n,:);
        
        % initialize new transition counts array 
        tr_counts_temp = zeros(nStates,nStates);
        
        % initialize state occupancy counts 
        state_counts_temp = zeros(1,nStates);
        
        % initialize array to store updated state counts
%         sc_slice = sample_chains_prob(:,n,:,:);
        
        for u = 1:us_factor_fine-1
            A_array_fw(:,:,u) = A_temp^iter_vec_fw(u);
            A_array_bk(:,:,u) = A_temp^iter_vec_bk(u);
        end
        for to = 1:nStates
            for from = 1:nStates
                % generate K x T array with state probabilities
                fwd_array = permute(A_array_fw(:,from,:),[1 3 2]);
                bkd_array = permute(A_array_bk(to,:,:),[2 3 1]);
                fwd_array = [zeros(nStates,1) fwd_array];    
                fwd_array(from,1) = 1;
                bkd_array = [bkd_array zeros(nStates,1)];    
                bkd_array(to,end) = 1;

                ss_slice = fwd_array(:,2:end).*bkd_array(:,1:end-1);%;
                ss_slice = ss_slice./ sum(ss_slice,1);            
                st_array = [zeros(nStates,1) ss_slice zeros(nStates,1)];    
                st_array(from,1) = 1;
                st_array(to,end) = 1;

                % calculate approximate fraction of time spent in each state
                state_counts = sum(st_array(:,1:end-1),2)/us_factor_fine;

                % calculate predicted state transition counts
                tr_counts_fw = repmat(A_temp,1,1,us_factor_fine).*permute(fwd_array,[3,1,2]);
                tr_counts = tr_counts_fw.*permute(bkd_array,[1,3,2]);
                tr_counts = sum(tr_counts./ sum(sum(tr_counts,1),2),3);
                
                % add estimates to arrays
                tr_counts(eye(nStates)==1) = 0;
                i_flags = from_array==from & to_array==to;
                
                n_tr = sum(i_flags(:));
                
                tr_counts_temp = tr_counts_temp + n_tr*tr_counts;
                state_counts_temp = state_counts_temp + n_tr*state_counts';
                
                for s = 1:nStates
                    sc_slice = sample_chains_prob{s}(1:end-1,n,:);
                    sc_slice(i_flags) = state_counts(s);
                    sample_chains_prob{s}(1:end-1,n,:) = sc_slice;
                    ec = sample_chains_prob{s}(end,n,:);
                    sample_chains_prob{s}(end,n,ec==s) = 1;
                end
            end
        end                  
        mcmcInfo.state_counts(:,:,n) = state_counts_temp;
        mcmcInfo.transition_count_array(:,:,n) = tr_counts_temp;
    end  
    mcmcInfo.sample_chains_prob = sample_chains_prob;
%     mcmcInfo.state_counts = sum(mcmcInfo.transition_count_array,1);