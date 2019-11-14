function logL_avg = calculate_likelihood(emissions_cell,pi_cell,v_curr)

% initialize likelihood vectors
logL_avg_vec = [];
logL_wt = [];
K = numel(v_curr);
T = numel(emissions_cell{1});
% iterate through traces
for e = 1:numel(pi_cell)
    e_vec = emissions_cell{e};
    pi_array = pi_cell{e};
    % get emissions provs for each observation/state
    emit_probs = poisspdf(repelem(e_vec,K),repmat(v_curr',1,T));
    % calculate avg state occupancy per timepoint        
    logL_avg_main = pi_array .* reshape(emit_probs,K,T);
%     logL_avg_start = log(pi0_curr(pt_array(:,1))) + log(v_curr(pt_array(:,1)));
    logL_avg_vec = [logL_avg_vec nanmean(logL_avg_main)];
    logL_wt = [logL_wt size(pt_array,2)];
end
logL_avg = sum(logL_avg_vec.*logL_wt) / sum(logL_wt);