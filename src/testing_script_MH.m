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
n_particles = 100;                 
n_steps = 100;
logL_vec = NaN(n_steps,1);
A_inf_array = NaN(K,K,n_steps);
v_inf_array = NaN(K,n_steps);
%% simple MCMC 
% specify hyperparameters
v_sigma = .2;
n_a = 100;

% initialize
A_curr = rand(K);
A_curr = A_curr ./ sum(A_curr);
v_curr = sort(rand(K,1))*2;
[V, D] = eig(A_curr);
[~, mi] = max(diag(D));
pi0_curr = V(:,mi)/sum(V(:,mi));
% calculate initial likelihood
[v_inf, A_inf, pi_cell] = particle_mcmc(n_particles,emissions_cell,A_curr,v_curr,pi0_curr);    
% calculate initial likelihood                        
logL_avg = calculate_likelihood(emissions_cell,pi_cell,v_curr);
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