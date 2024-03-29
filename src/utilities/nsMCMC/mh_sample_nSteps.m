function mcmcInfo = mh_sample_nSteps(mcmcInfo)

% pull out parameters for conveneice
% n_chains = mcmcInfo.n_chains;
% n_chains_eff = mcmcInfo.n_chains_eff;

% propose new nSteps values
%%
% tic
for iter = 1:mcmcInfo.nStep_tries_per_run
    nStepsCurr = mcmcInfo.nStepsCurr;
    ubVec = mcmcInfo.nStepsMax-nStepsCurr;
    lbVec = mcmcInfo.nStepsMin-nStepsCurr;
    nStepsProp = nStepsCurr + trandn(lbVec,ubVec)*mcmcInfo.nStepsPropSize;

    % generate MS2 kernels
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
end
% toc
%%
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