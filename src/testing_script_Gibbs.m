clear 
close all

addpath('utilities')

% define parameters
K =2;
R = [-.02, .04; .02 -.04];
tres = 3;
A = expm(R*tres);
v = [.1,2]';
seq_length = 40*60/3;
w = 7;
[V, D] = eig(A);
[~, mi] = max(diag(D));
pi0 = V(:,mi)/sum(V(:,mi));
noise = 1;
alpha = 1.4;
n_traces = 10;

emissions_cell = cell(1,n_traces);
for i = 1:n_traces
    synthetic_data = synthetic_prob_poisson(seq_length, alpha, K, w, A, ...
                                             v, noise, pi0);
    emissions_cell{i} = synthetic_data.loading_events;
end

%% Gibbs MCMC 

% basic inference params 
n_particles = 100; % number of particles to simulate
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
vp_mean = [.1, 1]; % prior on RNAP initiation rates
theta_prior = [2 2]./(n_particles); % 1 over number of intervals observed (make this of order N counts per round)
k_prior = vp_mean ./ theta_prior; % prior on k param


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_curr = sample_A_dirichlet(A_alpha, zeros(K));
v_curr = [gamrnd(k_prior(1),theta_prior(1),1) gamrnd(k_prior(2),theta_prior(2),1)];
% calculate pi0 
[V, D] = eig(A_curr);
[~, mi] = max(diag(D));
pi0_curr = V(:,mi)/sum(V(:,mi));
% calculate initial likelihood
tic
[v_cts_out, v_wts_out, A_count, pi_cell] = particle_mcmc(n_particles,emissions_cell,A_curr,v_curr,pi0_curr);    
toc
% calculate initial likelihood    
tic
[fwd_probs, logL_seq, logL_tot] = fwd_algorithm(emissions_cell, A_curr, v_curr, pi0_curr);
toc
% record
logL_vec(1) = logL_avg;
A_inf_array(:,:,1) = A_curr;
v_inf_array(:,1) = v_curr;
for s = 2:n_steps
    tic
    A_curr = A_inf_array(:,:,s-1);
    v_curr = v_inf_array(:,s-1);
    logL_curr = logL_vec(s-1);
    [V, D] = eig(A_curr);
    [~, mi] = max(diag(D));
    pi0_curr = V(:,mi)/sum(V(:,mi));
    % sequential MCMC step for sampling promoter trajectories   
    [v_inf, A_inf, pi_cell] = particle_mcmc(n_particles,emissions_cell,A_curr,v_curr,pi0_curr);                   
    % standard MCMC move for A
    A_prop = sample_A(A_curr,n_a);
    [V, D] = eig(A_prop);
    [~, mi] = max(diag(D));
    pi0_prop = V(:,mi)/sum(V(:,mi));
    logL_A_prop = calculate_likelihood(emissions_cell,pt_cell,v_curr);    
    % compare likelihoods
    take_A_move = logL_curr/logL_A_prop > rand();
    if take_A_move
        A_curr = A_prop;
        pi0_curr = pi0_prop;
        logL_curr = logL_A_prop;
    end
    % standard MCMC move for v
    v_prop = sample_v(v_curr,v_sigma);
    logL_v_prop = calculate_likelihood(emissions_cell,pt_cell,v_prop);    
    % compare likelihoods
    take_v_move = logL_curr/logL_v_prop > rand();
    if take_v_move
        v_curr = v_prop;        
        logL_curr = logL_v_prop;
    end    
    % record
    A_inf_array(:,:,s) = A_curr;
    v_inf_array(:,s) = v_curr;
    logL_vec(s) = logL_curr;
    toc
end