function mcmcInfo = initialize_chains(mcmcInfo)

% extract parameters
A = mcmcInfo.A_curr;
nStates = size(A,1);
pi0 = mcmcInfo.pi0_curr;
sql = mcmcInfo.seq_length;
n_chains = mcmcInfo.n_chains;
n_traces = mcmcInfo.n_traces;
state_options = 1:nStates;

% initialize
mcmcInfo.sample_chains = NaN(sql,n_chains,n_traces);

% draw initial state
for n = 1:n_traces
    mcmcInfo.sample_chains(1,:,n) = randsample(state_options,n_chains,true,pi0);    
    % simulate remaining time step
    for t = 2:sql
        % simulate transitions
        A_array = cumsum(A(:,mcmcInfo.sample_chains(t-1,:,n)));
        rand_array = repmat(rand(1,n_chains),nStates,1);
        mcmcInfo.sample_chains(t,:,n) = sum(rand_array > A_array) + 1;
    end    
end   
