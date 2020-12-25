clear 
close all

addpath('utilities')

% define parameters
R = [-.02, .04, 0; .02 -.05 .08; 0 .01 -.08];
tres = 3;
A = expm(R*tres);
K = size(A,1);
v = [.05,2,4]';
seq_length = 40*60/3;
nStates = 7;
[V, D] = eig(A);
[~, mi] = max(diag(D));
pi0 = V(:,mi)/sum(V(:,mi));
noise = 1;
alpha = 1.4;
n_traces = 10;
eps = 1e-2;

emissions_cell = cell(1,n_traces);
master_sim_struct = struct;
for i = 1:n_traces
    synthetic_data = synthetic_prob_poisson(seq_length, alpha, K, nStates, A, ...
                                             v, noise, pi0);
    emissions_cell{i} = synthetic_data.loading_events;
    
    fieldNames = fieldnames(synthetic_data);
    for f = 1:length(fieldNames)
        master_sim_struct(i).(fieldNames{f}) = synthetic_data.(fieldNames{f});
    end
end

%% Gibbs MCMC 

% basic inference params 
n_steps = 100; % number of MCMC steps (need to add convergence criteria)

% initialize arrays
logL_vec = NaN(n_steps,1);
A_inf_array = NaN(K,K,n_steps);
v_inf_array = NaN(K,n_steps);

% specify hyperparameters

% A prior--assume strongly diagonal PDF given short timescale
% take A columns to follow multinomial Dirichlet distribution
A_alpha = ones(K);%*n_particles*n_traces;
A_alpha(eye(K)==1) = A_alpha(eye(K)==1)*10; % distribution hyper params

% emission rate priors. For now assume v_i's to be Poisson. Thus, the
% conjugate prior is gamma distributed
vp_mean = [.1, 2.1 4.2]; % prior on RNAP initiation rates
theta_prior = [2 2 2]./(100); % 1 over number of intervals observed (make this of order N counts per round)
k_prior = vp_mean ./ theta_prior; % prior on k param


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_curr = sample_A_dirichlet(A_alpha, zeros(K));
for k = 1:K
    v_curr(k) = gamrnd(k_prior(k),theta_prior(k),1);
end

% calculate pi0 
[V, D] = eig(A_curr);
[~, mi] = max(diag(D));
pi0_curr = V(:,mi)/sum(V(:,mi));

% calculate forward and backward trajectories
[fwd_probs, logL_seq, logL_curr] = fwd_algorithm(emissions_cell, A_curr, v_curr, pi0_curr);

% record
logL_vec(1) = logL_curr;
A_inf_array(:,:,1) = A_curr;
v_inf_array(:,1) = v_curr;

step = 2;
deltaL = Inf;

while step < n_steps && deltaL > eps
  
    % get current parameter values
    A_curr = A_inf_array(:,:,step-1);
    v_curr = v_inf_array(:,step-1)';    
    [V, D] = eig(A_curr);
    [~, mi] = max(real(diag(D)));
    pi0_curr = V(:,mi)/sum(V(:,mi));
    
    % sequential MCMC step for sampling promoter trajectories  
    [fwd_probs, logL_seq, logL_curr] = fwd_algorithm(emissions_cell, A_curr, v_curr, pi0_curr);
    [bkd_probs] = bkd_algorithm(emissions_cell, A_curr, v_curr, pi0_curr);

    [n_init, n_steps, A_counts, posterior_prob_cell] = det_event_counts(...
        emissions_cell, fwd_probs, bkd_probs, A_curr, v_curr);
      
    % Draw A and v samples
    %%% A
    A_alpha_update = A_alpha + A_counts;
    A_curr = sample_A_dirichlet(A_alpha_update, zeros(K));
    
    %%% v
    theta_update = 1 ./ (theta_prior.^-1 + n_steps');
    k_update = k_prior + n_init';
    for k = 1:K
        v_curr(k) = gamrnd(k_update(k),theta_update(k),1);
    end    
    
    % record
    A_inf_array(:,:,step) = A_curr;
    v_inf_array(:,step) = v_curr;
    logL_vec(step) = logL_curr;    
    
    
    % calculate interation stats
    deltaL = Inf;
    if step > 5
        logLMean = movmean(logL_vec,5,'omitnan');
        logLDiffMean = movmean([0 diff(logL_vec)'],5,'omitnan');
        deltaL = logLDiffMean(step)/logLMean(step);
    end
    step = step + 1;
end