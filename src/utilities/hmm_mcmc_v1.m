function [A_inf_array, v_inf_array, logL_vec] = hmm_mcmc_v1(...
                    K,emissions_cell,n_steps,varargin)
                  
R = [-.02, .04, 0; .02 -.05 .08; 0 .01 -.08];
tres = 3;
A = expm(R*tres);
% set defaults
% A prior--assume strongly diagonal PDF given short timescale
% take A columns to follow multinomial Dirichlet distribution
A_prior = ones(K);%*n_particles*n_traces;
A_prior(eye(K)==1) = A_prior(eye(K)==1)*10; % distribution hyper params

% emission rate priors. For now assume v_i's to be Poisson. Thus, the
% conjugate prior is gamma distributed
vp_prior = [0.05 1.5 3]; % prior on RNAP initiation rates
target_n = 100;
for i = 1:numel(varargin)-1
    if ischar(varargin{i})
        eval([varargin{i} ' = varargin{i+1};']);
    end
end

% calculate weigth factor for gibbs sampling
n_obs = numel([emissions_cell{:}])/K;
a_factor = target_n / n_obs;

theta_prior = [1 1 1]/(numel([emissions_cell{:}])/K)/a_factor * 10;% * (a_factor*10))); % 1 over number of intervals observed (make this of order N counts per round)
k_prior = vp_prior ./ theta_prior; % prior on k param

% initialize arrays
logL_vec = NaN(n_steps,1);
A_inf_array = NaN(K,K,n_steps);
v_inf_array = NaN(K,n_steps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_curr = sample_A_dirichlet(A_prior, zeros(K));
for k = 1:K
    v_curr(k) = gamrnd(k_prior(k),theta_prior(k),1);
end
% calculate pi0 
[V, D] = eig(A_curr);
[~, mi] = max(diag(D));
pi0_curr = V(:,mi)/sum(V(:,mi));
% calculate forward and backward trajectories
[~, ~, logL_curr] = fwd_algorithm(emissions_cell, A_curr, v_curr, pi0_curr);

% record
logL_vec(1) = logL_curr;
A_inf_array(:,:,1) = A_curr;
v_inf_array(:,1) = v_curr;
for s = 2:n_steps    
%     tic
    % get current parameter values
    A_curr = A_inf_array(:,:,s-1);
    v_curr = v_inf_array(:,s-1)';    
    [V, D] = eig(A_curr);
    [~, mi] = max(diag(D));
    pi0_curr = V(:,mi)/sum(V(:,mi));
    
    % sequential MCMC step for sampling promoter trajectories  
    [fwd_probs, ~, logL_curr] = fwd_algorithm(emissions_cell, A_curr, v_curr, pi0_curr);
    [bkd_probs] = bkd_algorithm(emissions_cell, A_curr, v_curr, pi0_curr);

    [n_init, n_steps, A_counts, ~] = det_event_counts(...
        emissions_cell, fwd_probs, bkd_probs, A_curr, v_curr);
    % Draw A and v samples
    %%% A
    A_alpha_update = A_counts + A_prior;
    A_curr = sample_A_dirichlet(A_alpha_update, zeros(K));
    %%% v
    theta_update = 1 ./ ( n_steps');
    k_update = k_prior + n_init';
    for k = 1:K
        v_curr(k) = gamrnd(k_update(k),theta_update(k),1);
    end    
    % record
    A_inf_array(:,:,s) = A_curr;%A_counts ./ sum(A_counts);%A_curr;
    v_inf_array(:,s) = v_curr;%n_init ./ n_steps;%v_curr;
    logL_vec(s) = logL_curr;    
%     toc
end