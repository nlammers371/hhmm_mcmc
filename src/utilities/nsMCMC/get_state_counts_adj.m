function [mcmcInfo, F_array, y_array] = get_state_counts_adj(mcmcInfo)

seq_length = mcmcInfo.seq_length;
n_traces = mcmcInfo.n_traces;
nStates = mcmcInfo.nStates;
n_chains = mcmcInfo.n_chains;
coeff_MS2_us = mcmcInfo.coeff_MS2_us; 
us_factor = mcmcInfo.upsample_factor; 

% generate F count arrays
F_array = zeros(seq_length*n_traces,nStates,n_chains);        
y_array = NaN(seq_length*n_traces,n_chains);    

for c = 1:n_chains
    for n = 1:n_traces        
        ind1 = (n-1)*seq_length+1;
        ind2 = n*seq_length;

        % record observed fluo                        
        y_array(ind1:ind2,c) = mcmcInfo.observed_fluo(:,n);

        for m = 1:nStates
            % record counts
            state_counts = convn(coeff_MS2_us(:,c),mcmcInfo.sample_chains_prob{m}(:,c,n),'full')/us_factor;            
            state_counts = state_counts(1:end-size(coeff_MS2_us,1)+1,:);
            state_counts = state_counts(us_factor:us_factor:end,:);
            F_array(ind1:ind2,m,c) = state_counts;                        
        end
    end
end 