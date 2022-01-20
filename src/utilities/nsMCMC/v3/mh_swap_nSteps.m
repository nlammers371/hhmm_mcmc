function mcmcInfo = mh_swap_nSteps(mcmcInfo)


% generate distribution from previous N steps (excluding burn-in)
start_i = ceil((mcmcInfo.step-mcmcInfo.swapInc+1)/mcmcInfo.update_increment);
stop_i = floor(mcmcInfo.step/mcmcInfo.update_increment);

nStepsDist = mcmcInfo.n_steps_inf_array(start_i:stop_i,:);
nStepsDist = nStepsDist(:);
nStepsWeight = exp(mcmcInfo.logL_vec(start_i:stop_i,:));
nStepsWeight = nStepsWeight(:);
nStepsProp = randsample(nStepsDist,mcmcInfo.n_chains,true,nStepsWeight);

% randomly sample MS2 kernels
coeff_MS2_prop = NaN(size(mcmcInfo.coeff_MS2));
for n = 1:length(nStepsProp)
    alpha = mcmcInfo.alpha_frac*nStepsProp(n);
    coeff_MS2_prop(:,n) = ms2_loading_coeff_frac(alpha, nStepsProp(n), mcmcInfo.nStepsMax);
end

[mcmcInfo, trace_logL_array, trace_temp_array] = calculateNStepLogL(mcmcInfo, coeff_MS2_prop);

% calculate MH metric
L_factor = exp(sum(trace_logL_array));

% perform MH move
mh_flags = L_factor >= rand(size(L_factor));

% update    
mcmcInfo.coeff_MS2(:,mh_flags) = coeff_MS2_prop(:,mh_flags);
mcmcInfo.nStepsCurr(mh_flags) = nStepsProp(mh_flags);
mcmcInfo.sample_fluo(:,mh_flags,:) = trace_temp_array(:,mh_flags,:);

% if mcmcInfo.temperingFlag
%     
%     for n = 1:mcmcInfo.n_chains
%         % propose memory swaps if we're doing tempering
%         L_factor = exp((total_log_likelihoods(1,:,3)-total_log_likelihoods(1,:,1))./temperature_array(1,:) + ...
%                                 (total_log_likelihoods(1,:,4)-total_log_likelihoods(1,:,2))./temperature_array(2,:));
% 
%         % perform MH move
%         mh_flags = L_factor >= rand(size(L_factor));
%     end
% end                      


% store result if appropriate
if mcmcInfo.update_flag
    update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    mcmcInfo.n_steps_inf_array(update_index,:) = mcmcInfo.nStepsCurr;
end