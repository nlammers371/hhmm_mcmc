function to_array = allowed_to_array(state, K, w)   
    % List of all compound states that can transition into the given 
    % array of compound states.
    % 
    % INPUTS
    % state: compound state array
    % K: number of naive states
    % w: memory
    % 
    % OUTPUTS
    % l : the list of all compound states that can transition into the 
    %     given compound state. Note: the fact that should be an overlap of 
    %     (w-1) naive states in the past and present is used in the 
    %     calculations.
    
    first_naive = floor( (state-1) ./ K^(w-1));
    tail_naive = state-1 - first_naive * K^(w-1);
%     l = 1 + (K*tail_naive : (K*tail_naive + K - 1));
    
    % generate new array
    sz = size(state);
    sz(1) = K;
    to_array = NaN(sz);
    for k = 1:K
        to_array(k,:) = K*tail_naive + k;
    end