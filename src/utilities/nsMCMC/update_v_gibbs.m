function mcmcInfo = update_v_gibbs(mcmcInfo, F_array, y_array)

n_chains = mcmcInfo.n_chains; 

if mcmcInfo.update_increment~=1
    update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
else
    update_index = mcmcInfo.step;
end

for c = 1:n_chains   
        
    M = ((F_array(:,:,c)'*F_array(:,:,c))) + 1;    
    b = ((F_array(:,:,c)'*y_array(:,c)));

    % calculate mean and variance
    v_lsq = M\b;
    if ~any(isnan(v_lsq))
        v_mean = (M + mcmcInfo.M0)^-1 * (mcmcInfo.M0*mcmcInfo.v0(c,:)' + M*v_lsq);
        v_cov_mat = inv(mcmcInfo.sigma_curr(c)^-2 * M +  mcmcInfo.sigma_curr(c)^-2 * inv(mcmcInfo.M0));
    else % if issue with lsq solution, draw from prior
         v_cov_mat = mcmcInfo.sigma_curr(c)^2 * inv(mcmcInfo.M0);
         v_mean = mcmcInfo.v0(c,:);
    end
    % sample            
    try
%             [mcmcInfo.v_curr(c,:), swap_indices(c,:)] = sort(mvnrnd(v_mean, v_cov_mat)'); 
        mcmcInfo.v_curr(c,:) = mvnrnd(v_mean, v_cov_mat)'; 
    catch
        error('check')
    end

    if mcmcInfo.update_flag
        mcmcInfo.v_inf_array(update_index,:,c) = mcmcInfo.v_curr(c,:);   
    end
end