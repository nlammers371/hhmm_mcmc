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
w = 7*6;
[V, D] = eig(A);
[~, mi] = max(diag(D));
pi0 = V(:,mi)/sum(V(:,mi));
noise = 1;
alpha = 1.4;
n_traces = 20;
n_steps = 50;

emissions_cell = cell(1,n_traces);
fluo_cell = cell(1,n_traces);
for i = 1:n_traces
    synthetic_data = synthetic_prob_poisson(seq_length, alpha, K, w, A, ...
                                             v, noise, pi0);
    emissions_cell{i} = synthetic_data.loading_events;
    fluo_cell{i} = synthetic_data.fluo_MS2;
end

%%
f_max = max([fluo_cell{:}]);
n_step_approx = w;
v_max_guess = f_max / n_step_approx * 1.3
vp_prior = [.05, v_max_guess/2, v_max_guess]

%%
n_runs = 25;
burn_in = 10;
A_inf_cell = cell(1,n_runs);
v_inf_cell = cell(1,n_runs);
logL_cell = cell(1,n_runs);
logL_vec = NaN(1,n_runs);

disp('conducting mcmc inference...')
for n = 1:n_runs
    vp_prior = sort(rand(1,3)*v_max_guess);    
    tic
    [A_inf_cell{n}, v_inf_cell{n}, logL_cell{n}] = hmm_mcmc_v1(...
                        K,emissions_cell,n_steps,'vp_prior',vp_prior);
    disp(['done (' sprintf('%02d',n) ' of ' num2str(n_runs) ')'])
    logL_vec(n) = mean(logL_cell{n}(end-burn_in+1:end));
    toc
end