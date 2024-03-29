function mcmcInfo = adjust_tr_counts(mcmcInfo)

    nStates = mcmcInfo.nStates;
    Q_curr = mcmcInfo.Q_curr;
    n_chains = mcmcInfo.n_chains;
    r_mat = rand(size(Q_curr))/1e5;
    d_mat = repmat(eye(nStates),1,1,n_chains)==1;
    r_mat(d_mat==1) = 0;
    r_mat(d_mat==1) = -sum(r_mat,1);
    
    Q_curr = Q_curr + r_mat; % to avoid divide by 0 errors
    
    
    
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
    elseif nStates == 3        
        % extract rates
        k21 = reshape(Q_curr(2,1,:),1,[]);
        k32 = reshape(Q_curr(3,2,:),1,[]);
        k23 = reshape(Q_curr(2,3,:),1,[]);
        k12 = reshape(Q_curr(1,2,:),1,[]);
        
        % extract initial counts
        counts_raw = mcmcInfo.transition_count_array;
        
        % time between samples
        T = mcmcInfo.tres/mcmcInfo.upsample_factor;        
        
        % calculate path weights
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Column 1 of A
        
        % 1->1
        w011 = 1-calculate_jump_weights1(k21,T); %1->1
        w211 = calculate_jump_weights2(k21,k12,T); % 1->2->1
        w411Rep = zeros(size(w211));%
        w411 = zeros(size(w211));
        if mcmcInfo.rateSamplingHRFlag
            w411Rep = calculate_jump_weights4Rep(k21,k12,T); % 1->2->1->2->1
            w411 = calculate_jump_weights4(k21,k32,k23,k12,T); % 1->2->3->2->1
        end
        denom11 = w011+w211+w411+w411Rep;
        
        % 1->2        
        w121 = calculate_jump_weights1(k21,T); %1->2
        w321Rep = zeros(size(w211));%calculate_jump_weights3Rep(k21,k12,T); % 1->2->1->2
        w321 = zeros(size(w211));%calculate_jump_weights3(k21,k32,k23,T); % 1->2->3->2
        if mcmcInfo.rateSamplingHRFlag
            w321Rep = calculate_jump_weights3Rep(k21,k12,T); % 1->2->1->2
            w321 = calculate_jump_weights3(k21,k32,k23,T); % 1->2->3->2
        end
        denom21 = w121+w321+w321Rep;
        
        % 1->3 
        w231 = calculate_jump_weights2(k21,k32,T); % 1->2->3        
        w431_1 = zeros(size(w231));%calculate_jump_weights4Rep1(k21,k12,k32,T);%calculate_jump_weights2(k21,k32,T); % 1->2->1->2->3        
        w431_2 = zeros(size(w231));%calculate_jump_weights4Rep2(k21,k32,k23,T);%calculate_jump_weights2(k21,k32,T); % 1->2->3->2->3
        if mcmcInfo.rateSamplingHRFlag
            w431_1 = calculate_jump_weights4Rep1(k21,k12,k32,T);%calculate_jump_weights2(k21,k32,T); % 1->2->1->2->3        
            w431_2 = calculate_jump_weights4Rep2(k21,k32,k23,T);%calculate_jump_weights2(k21,k32,T); % 1->2->3->2->3
        end
        denom31 = w231+w431_1+w431_2;
        
        % 2->2
        w022 = 1-calculate_jump_weights1(k12+k32,T); %2->2
        w222_1 = calculate_jump_weights2(k12,k21,T); % 2->1->2
        w222_2 = calculate_jump_weights2(k32,k23,T); % 2->3->2
        w422_1 = zeros(size(w231));%calculate_jump_weights4Rep(k12,k21,T); % 2->1->2->1->2
        w422_2 = zeros(size(w231));%calculate_jump_weights4Rep(k32,k23,T); % 2->3->2->3->2
        w422_3 = zeros(size(w231));%calculate_jump_weights4(k12,k21,k32,k23,T); % 2->1->2->3->2
        w422_4 = zeros(size(w231));%calculate_jump_weights4(k32,k23,k12,k21,T); % 2->3->2->1->2
        if mcmcInfo.rateSamplingHRFlag
            w422_1 = calculate_jump_weights4Rep(k12,k21,T); % 2->1->2->1->2
            w422_2 = calculate_jump_weights4Rep(k32,k23,T); % 2->3->2->3->2
            w422_3 = calculate_jump_weights4(k12,k21,k32,k23,T); % 2->1->2->3->2
            w422_4 = calculate_jump_weights4(k32,k23,k12,k21,T); % 2->3->2->1->2
        end
        denom22 = w022+w222_1+w222_2+w422_1+w422_2+w422_3+w422_4;
        
        % 2->1
        w112 = calculate_jump_weights1(k12,T); %2->1
        w312_1 = zeros(size(w211));%calculate_jump_weights3(k32,k23,k12,T); % 2->3->2->1
        w312_2 = zeros(size(w211));%calculate_jump_weights3Rep(k12,k21,T); % 2->1->2->1
        if mcmcInfo.rateSamplingHRFlag
            w312_1 = calculate_jump_weights3(k32,k23,k12,T); % 2->3->2->1
            w312_2 = calculate_jump_weights3Rep(k12,k21,T); % 2->1->2->1
        end
        denom12 = w112+w312_1+w312_2;
        
        % 2->3
        w132 = calculate_jump_weights1(k32,T); %2->3
        w332_1 = zeros(size(w211));%calculate_jump_weights3Rep(k32,k23,T); % 2->3->2->3
        w332_2 = zeros(size(w211));%calculate_jump_weights3(k12,k21,k32,T); % 2->1->2->3
        if mcmcInfo.rateSamplingHRFlag
            w332_1 = calculate_jump_weights3Rep(k32,k23,T); % 2->3->2->3
            w332_2 = calculate_jump_weights3(k12,k21,k32,T); % 2->1->2->3
        end
        denom32 = w132+w332_1+w332_2;
        
        % 3->3
        w033 = 1-calculate_jump_weights1(k23,T); %3->3
        w233 = calculate_jump_weights2(k23,k32,T); % 3->2->3       
        w433_1 = zeros(size(w231));%calculate_jump_weights4Rep(k23,k32,T); % 3->2->3->2->3
        w433_2 = zeros(size(w231));%calculate_jump_weights4(k23,k12,k21,k32,T); % 3->2->1->2->3
        if mcmcInfo.rateSamplingHRFlag
            w433_1 = calculate_jump_weights4Rep(k23,k32,T); % 3->2->3->2->3
            w433_2 = calculate_jump_weights4(k23,k12,k21,k32,T); % 3->2->1->2->3
        end
        denom33 = w033+w233+w433_1+w433_2;
        
        % 3->1
        w213 = calculate_jump_weights2(k23,k12,T); % 3->2->1        
        w413_1 = zeros(size(w231));%calculate_jump_weights4Rep1(k23,k32,k12,T);%calculate_jump_weights2(k21,k32,T); % 3->2->3->2->1        
        w413_2 = zeros(size(w231));%calculate_jump_weights4Rep2(k23,k21,k12,T);%calculate_jump_weights2(k21,k32,T); % 3->2->1->2->1
        if mcmcInfo.rateSamplingHRFlag
            w413_1 = calculate_jump_weights4Rep1(k23,k32,k12,T);%calculate_jump_weights2(k21,k32,T); % 3->2->3->2->1        
            w413_2 = calculate_jump_weights4Rep2(k23,k21,k12,T);%calculate_jump_weights2(k21,k32,T); % 3->2->1->2->1
        end
        denom13 = w213+w413_1+w413_2;
        
        % 3->2
        w123 = calculate_jump_weights1(k23,T); %3->2
        w323_1 = zeros(size(w211));%calculate_jump_weights3Rep(k23,k32,T); % 3->2->3->2
        w323_2 = zeros(size(w211));%calculate_jump_weights3(k23,k12,k21,T); % 3->2->1->2
        if mcmcInfo.rateSamplingHRFlag
            w323_1 = calculate_jump_weights3Rep(k23,k32,T); % 3->2->3->2
            w323_2 = calculate_jump_weights3(k23,k12,k21,T); % 3->2->1->2
        end
        denom23 = w123+w323_1+w323_2;
        
        % apply adjustments        
        counts_adjusted = zeros(size(counts_raw));
        
        counts_raw_rs = permute(reshape(counts_raw,nStates^2,1,n_chains),[1 3 2]);
        
        % generate weight array (easier for bookkeeping)
        weight_array = [[w011 ; w211 ; w411Rep ;  w411]./denom11;...
                               [ w121 ;  w321 ; w321Rep]./denom21;...
                               [ w231 ;  w431_1 ; w431_2]./denom31;...
                               [ w112 ;  w312_1 ; w312_2]./denom12;...
                               [ w022 ; w222_1 ; w222_2 ; w422_1 ; w422_2 ; w422_3 ; w422_4]./denom22;...
                               [ w132 ; w332_1 ; w332_2]./denom32;...
                               [ w213 ; w413_1 ; w413_2]./denom13;...
                               [ w123 ; w323_1 ; w323_2]./denom23;...
                               [ w033 ; w233 ; w433_1 ; w433_2]./denom33];
        
        lin_ind_vec = [1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 5 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 9]';
        
        % generate array of counts
                        % 1->2,2->1,3->2,2->3
        tr_count_array =   [0   0    0   0 ; % w011: 1->1
                            1   1    0   0 ;% w211: 1->2->1
                            2   2    0   0 ;% w411Rep: 1->2->1->2->1
                            1   1    1   1 ;% w411: 1->2->3->2->1
                            1   0    0   0 ;% w121: 1->2
                            1   0    1   1 ;% w321_1: 1->2->3->2
                            2   1    0   0 ;% w321Rep: 1->2->1->2
                            1   0    0   1 ;% w231: 1->2->3
                            2   1    0   1 ;% w431_1: 1->2->1->2->3  
                            1   0    1   2 ;% w431_2: 1->2->3->2->3  
                            0   1    0   0 ;% w112: 2->1
                            0   1    1   1 ;% w312_1: 2->3->2->1
                            1   2    0   0 ;% w312_2: 2->1->2->1
                            0   0    0   0 ;% w022: 
                            1   1    0   0 ;%w222_1: 2->1->2
                            0   0    1   1 ;%w222_2: 2->3->2
                            2   2    0   0 ;%w422_1: 2->1->2->1->2
                            0   0    2   2 ;%w422_2: 2->3->2->3->2
                            1   1    1   1 ;%w422_3: 2->1->2->3->2
                            1   1    1   1 ;%w422_4: 2->3->2->1->2
                            0   0    0   1 ;%w132: 2->3
                            0   0    1   2 ;%w332_1: 2->3->2->3
                            1   1    0   1 ;%w332_2: 2->1->2->3
                            0   1    1   0 ;%w213: 3->2->1        
                            0   1    2   1 ;%w413_1: 3->2->3->2->1        
                            1   2    1   0 ;%w413_2: 3->2->1->2->1
                            0   0    1   0 ;%w123: 3->2
                            0   0    2   1 ;%w323_1: 3->2->3->2
                            1   1    1   0 ;%w323_2: 3->2->1->2
                            0   0    0   0 ;%w033: 3->3
                            0   0    1   1 ;%w233: 3->2->3       
                            0   0    2   2 ;%w433_1: 3->2->3->2->3
                            1   1    1   1 ];%w433_2: 3->2->1->2->3
                     
                          
        % perform re-weighting calcultions
        weight_array_ct = weight_array.*counts_raw_rs(lin_ind_vec,:);
                
        counts_adjusted(2,1,:) = sum(weight_array_ct.*tr_count_array(:,1),1);
        counts_adjusted(1,2,:) = sum(weight_array_ct.*tr_count_array(:,2),1);
        counts_adjusted(2,3,:) = sum(weight_array_ct.*tr_count_array(:,3),1);
        counts_adjusted(3,2,:) = sum(weight_array_ct.*tr_count_array(:,4),1);        
        
        mcmcInfo.transition_count_array = counts_adjusted;
%         toc        
    end
    
    