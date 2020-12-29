function logL_fluo = calculate_fluo_logL_full(mcmcInfo)
                  
    % extrac useful parameters 
    sigma_curr = mcmcInfo.sigma_curr;    
    coeff_MS2 = mcmcInfo.coeff_MS2;                
    
    % extract relevant promoter state fragment
    initiation_array = mcmcInfo.v_curr(mcmcInfo.sample_chains); 
                                        
    % calculate predicted fluorescence
    fluo_array = convn(coeff_MS2,initiation_array,'full');             
    fluo_array = fluo_array(1:end-length(coeff_MS2)+1,:,:);%+1(ind-prevInd+1:end-nSteps+1,:,:);
    
    % get differences
    logL_fluo_full = ((permute(mcmcInfo.observed_fluo,[1 3 2])-fluo_array)./sigma_curr).^2;
                                      
    % take average
    logL_fluo = -mean(logL_fluo_full(:));%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);        
   