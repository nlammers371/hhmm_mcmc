% Script to experiment with non-sequential MCMC sampling for parameter
% inference
clear 
close all force

addpath(genpath('utilities'))

% initialize info structure
sampling_res = 16.85;
trueParams = setParamsReduced3state(sampling_res);

%%%%%%%%%%%%%%%%%%%%% Simulated data %%%%%%%%%%%%%%%%
% basic inference params 
n_chains = 50;
n_traces = 20;
seq_length = 120;
inferMemory = 0;

% global mcmcInfo
mcmcInfo.n_mcmc_steps = 50; 
mcmcInfo.burn_in = 50;
mcmcInfo.resampleTracesFlag = 1;
mcmcInfo.rs_freq = 1;
mcmcInfo.tres = sampling_res;
mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run

% characteristics of simulated data
mcmcInfo.upsample_factor = 1;
mcmcInfo.n_traces = n_traces;
mcmcInfo.seq_length = seq_length; % length of simulated traces in time steps

mcmcInfo.mhInferenceFlag = 1;
mcmcInfo.reducedModelFlag = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = 7;
trueParams.nSteps = nSteps;
trueParams.discrete_data_flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trueParams.n_traces = n_traces;
trueParams.seq_length = seq_length;
trueParams = generateSimulatedData(trueParams);

%%%%%%%%%%%%%%%%%%%%% Initialize sampling struct %%%%%%%%%%%%%%%%

% add "known" hyperparameters
mcmcInfo.nStates = trueParams.nStates;    
mcmcInfo.alpha_frac = trueParams.alpha_frac;
mcmcInfo.observed_fluo = trueParams.observed_fluo;

if ~inferMemory
    mcmcInfo.nSteps = trueParams.nSteps;
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
mcmcInfo = setMCMCOptions(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize inference arrays and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initializeInferenceArrays(mcmcInfo);
mcmcInfo = initializeVariablesBasicRandom(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmcInfo = initialize_chains_v3(mcmcInfo);

% get predicted fluorescence
mcmcInfo = predict_fluo_full_v3(mcmcInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct full inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
mcmcInfo.total_loss_vec = NaN(mcmcInfo.n_mcmc_steps,1);
mcmcInfo.mh_move_vec = false(mcmcInfo.n_mcmc_steps-1,1);

max_k = 1/sampling_res;
max_r = 5/sampling_res;

true_params = [trueParams.kon trueParams.koff  trueParams.k_corr_factor ...
               trueParams.r trueParams.r_corr_factor trueParams.sigma];
             
% define prop sizes and boundaries
%             [kon koff k_coop  r r_coop   sigma]
param_lb_vec = [0   0  log(2/3) 0 log(2/3)  1];
param_ub_vec = [max_k max_k log(1.5) max_r log(1.5)  10];
prop_size_vec = 0.025*param_ub_vec;

% calculate initial logL
mcmcInfo.step = 1;

mcmcInfo = resample_chains_PF_LL(mcmcInfo);

mcmcInfo.total_loss_vec(1) = mcmcInfo.LL;%*numel(mcmcInfo.observed_fluo);

wb = waitbar(1/mcmcInfo.n_mcmc_steps,'Estimating parameter values using MH MCMC');

for step = 2:mcmcInfo.n_mcmc_steps
  
    waitbar(step/mcmcInfo.n_mcmc_steps,wb);
    
    mcmcInfo.step = step;
    
    % extract current parameters
    k_curr = mcmcInfo.k_curr;
    r_curr = mcmcInfo.r_curr;
    logL_curr = mcmcInfo.total_loss_vec(step-1);
    sample_fluo_curr = mcmcInfo.sample_fluo;
    sample_chains_curr = mcmcInfo.sample_chains;
    logL_vec_curr = mcmcInfo.logL_vec(step-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Propose new MH move
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    current_params = [k_curr r_curr mcmcInfo.sigma_curr(1)];
    lb_array = (param_lb_vec - current_params)./prop_size_vec;
    ub_array = (param_ub_vec - current_params)./prop_size_vec;
    new_params = current_params+prop_size_vec.*trandn(lb_array,ub_array)';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% calculate new likelihood
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % A
    R_Matrix = Q_helper_fun(new_params(1),new_params(2),new_params(3));
    mcmcInfo.A_curr = repmat(expm(R_Matrix*mcmcInfo.tres),1,1,n_chains);
    
    % v
    r_vec = r_helper_fun(new_params(4), new_params(5));
    mcmcInfo.v_curr = repmat(r_vec*mcmcInfo.tres,n_chains,1);
    
    % calculate pi0 
    [V, D] = eig(mcmcInfo.A_curr(:,:,1));
    [~, mi] = max(real(diag(D)));
    mcmcInfo.pi0_curr = repmat((V(:,mi)/sum(V(:,mi)))',n_chains,1);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Update latent states based upon current parameter estimates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     mcmcInfo = initialize_chains_v3(mcmcInfo);
% 
%     % get predicted fluorescence
%     mcmcInfo = predict_fluo_full_v3(mcmcInfo);

    mcmcInfo = resample_chains_PF_LL(mcmcInfo);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate current model likelihood
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    % calculate updated logL        
    logL_prop = mcmcInfo.LL;%*numel(mcmcInfo.observed_fluo);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% perform MH move
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    mh_val = exp(logL_prop-logL_curr);
    
    if mh_val > rand() % accept
        % update working variable values 
        mcmcInfo.k_curr = new_params(1:3);
        mcmcInfo.r_curr = new_params(4:5);
        mcmcInfo.sigma_curr = new_params(6);
        
        % update parameters
        mcmcInfo.k_inf_array(step,:) = new_params(1:3);
        mcmcInfo.r_inf_array(step,:) = new_params(4:5);
        mcmcInfo.sigma_inf_array(step) = new_params(6);
        
        % update total loss vec
        mcmcInfo.total_loss_vec(step) = logL_prop;
        
        mcmcInfo.mh_move_vec(step-1) = true;
        
    else % reject
      
        % keep previous values
        mcmcInfo.k_inf_array(step,:) = mcmcInfo.k_inf_array(step-1,:);
        mcmcInfo.r_inf_array(step,:) = mcmcInfo.r_inf_array(step-1,:);
        mcmcInfo.sigma_inf_array(step) = mcmcInfo.sigma_inf_array(step-1); 
        
        % update total loss vec
        mcmcInfo.total_loss_vec(step) = logL_curr;
        
        % revert logL array
        mcmcInfo.logL_vec(step,:) = logL_vec_curr;
        
        % revert fluorescence values
        mcmcInfo.sample_fluo = sample_fluo_curr;
        mcmcInfo.sample_chains = sample_chains_curr;
%         mcmcInfo.logL_vec(step,:) = mcmcInfo.logL_vec(step-1,:);  
    end    

end
delete(wb)
% loss_init = calculateModelLoss(params_init);
% loss_true = calculateModelLoss(params_true);
% % mcmcInfo = inferenceWrapper(mcmcInfo);  
% loss_fun = @(params) -calculateModelLoss(params);
% tic
% % options = optimoptions('fmincon','ConstraintTolerance',1);
% [params, loss,~,output] = fmincon(loss_fun,params_init,[],[],[],[],[1e-3 1e-3 -0.7 0 -0.7 0],[.1 0.1 0.7 5/mcmcInfo.tres 0.7 5],[]);%,options);
% toc