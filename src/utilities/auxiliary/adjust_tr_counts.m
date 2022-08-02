function mcmcInfo = adjust_tr_counts(mcmcInfo)

    nStates = mcmcInfo.nStates;
    Q_curr = mcmcInfo.Q_curr;
%     n_chains = mcmcInfo.n_chains;
    
    if nStates == 2
%         if mcmcInfo.step > 51
%             error('check')
%         end
%         tic
        k1 = reshape(Q_curr(2,1,:),1,[]);
        k2 = reshape(Q_curr(1,2,:),1,[]);
        counts_raw = mcmcInfo.transition_count_array;
        T = mcmcInfo.tres/mcmcInfo.upsample_factor;        
        % calculate path weights
        
        % 1->1
        w011 = 1-calculate_jump_weights1(k1,T); %1->1
        w211 = calculate_jump_weights2(k1,k2,T); % 1->2->1
        w411 = calculate_jump_weights4Rep(k1,k2,T); % 1->2->1->2->1
        denom11 = w011+w211+w411;
        
        % 1->2
        
        w121 = calculate_jump_weights1(k1,T); %1->2
        w321 = calculate_jump_weights3Rep(k1,k2,T); % 1->2->1->2
        denom21 = w121+w321;
        
        % 2->2
        w022 = 1-calculate_jump_weights1(k2,T); %2->2
        w222 = calculate_jump_weights2(k2,k1,T); % 2->1->2
        w422 = calculate_jump_weights4Rep(k2,k1,T); % 2->1->2->1->2
        denom22 = w022+w222+w422;
        
        % 2->1
        w112 = calculate_jump_weights1(k2,T); %2->1
        w312 = calculate_jump_weights3Rep(k2,k1,T); % 2->1->2->1
        denom12 = w112+w312;
        
        % apply adjustments        
        counts_adjusted = zeros(size(counts_raw));
        ct11 = reshape(counts_raw(1,1,:),1,[]);       
        ct21 = reshape(counts_raw(2,1,:),1,[]);       
        ct22 = reshape(counts_raw(2,2,:),1,[]);       
        ct12 = reshape(counts_raw(1,2,:),1,[]);      
        
        % 1->1         
        counts_adjusted(1,1,:) = (w011./denom11).*ct11;
        
        % 1-> 2        
        counts_adjusted(2,1,:) = (w211./denom11).*ct11 + 2.*(w411./denom11).*ct11 + (w121./denom21).*ct21 + 2.*(w321./denom21).*ct21 + ...
                                 (w222./denom22).*ct22 + 2.*(w422./denom22).*ct22 + (w312./denom12).*ct12;
        % 2->2         
        counts_adjusted(2,2,:) = (w022./denom22).*ct22;
        
        % 1-> 2        
        counts_adjusted(1,2,:) = (w211./denom11).*ct11 + 2.*(w411./denom11).*ct11 + (w112./denom12).*ct12 + 2.*(w312./denom12).*ct12 + ...
                                 (w222./denom22).*ct22 + 2.*(w422./denom22).*ct22 + (w321./denom21).*ct21;
        
        mcmcInfo.transition_count_array = counts_adjusted;
%         toc        
    end