function mcmcInfo = mh_update_A(mcmcInfo)

n_chains = mcmcInfo.n_chains;
nStates = mcmcInfo.nStates;
us_factor = mcmcInfo.upsample_factor; 
update_flag = mcmcInfo.update_flag;

if mcmcInfo.update_increment~=1
    update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
else
    update_index = mcmcInfo.step;
end

pi0_counts = mcmcInfo.state_counts/us_factor;

if mcmcInfo.nStates == 3
    lin_indices_to = repmat([2 4 6 8],n_chains,1) + (nStates^2)*(0:n_chains-1)';    
elseif mcmcInfo.nStates == 2
    lin_indices_to = repmat([2 3],n_chains,1) + (nStates^2)*(0:n_chains-1)';
end

% generate linear indices encoding "observed" state transitions
from_array = mcmcInfo.sample_chains(1:end-1,:,:);
to_array = mcmcInfo.sample_chains(2:end,:,:);        
row_col_array = (from_array-1)*nStates+ to_array + nStates^2*mcmcInfo.chain_id_ref;
lin_index_array = row_col_array;

% tic
for iter = 1:mcmcInfo.Q_n_tries_per_run
    % extract current rate matrix
    Q_curr = mcmcInfo.Q_curr;        
    % reshape for sampling
    Q_curr_array = Q_curr(lin_indices_to);
    % draw random proposals
    ubArray = mcmcInfo.QMax-Q_curr_array;
    lbArray = 0-Q_curr_array;
    prop_deltas = reshape(trandn(lbArray,ubArray)*mcmcInfo.QPropSize,n_chains,[]);
    % generate proposed rate matrix
    Q_prop = zeros(size(Q_curr));
    Q_prop(lin_indices_to) = Q_curr_array + prop_deltas;        
    Q_prop(repmat(eye(nStates)==1,1,1,n_chains)) = -sum(Q_prop);
    % convert to transition probability
    A_prop = NaN(size(Q_prop));
    for n = 1:n_chains
        A_prop(:,:,n) = expm(Q_prop(:,:,n)*mcmcInfo.tres/mcmcInfo.upsample_factor);
    end
    A_curr_log = log(mcmcInfo.A_curr);
    A_prop_log = log(A_prop);
    
    % calculate likelihoods
    logL_transition_curr = sum(sum(A_curr_log(lin_index_array),1),3);
    logL_transition_prop = sum(sum(A_prop_log(lin_index_array),1),3);    

    % calculate MH metric
    L_factor = exp(logL_transition_prop-logL_transition_curr);

    % perform MH move
    mh_flags = L_factor >= rand(size(L_factor));

    % update    
    for n = 1:n_chains
        if mh_flags(n)
            mcmcInfo.Q_curr(:,:,n) = Q_prop(:,:,n);
            mcmcInfo.A_curr(:,:,n) = A_prop(:,:,n);
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