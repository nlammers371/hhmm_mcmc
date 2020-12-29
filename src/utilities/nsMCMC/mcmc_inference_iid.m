function mcmcInfoIID = mcmc_inference_iid(mcmcInfoIID)
    disp('Obtaining initial parameter estimates...');
    mcmcInfoIID.iid_flag = 1;        
    mcmcInfoIID.A_curr = ones(mcmcInfoIID.nStates)./mcmcInfoIID.nStates;
    mcmcInfoIID.pi0_curr = ones(mcmcInfoIID.nStates,1)./mcmcInfoIID.nStates;
    mcmcInfoIID = initialize_chains(mcmcInfoIID);
    
    for step = 1:mcmcInfoIID.n_mcmc_steps %mcmcInfo.n_mcmc_steps                
        mcmcInfoIID.step = step;    
        % resample chains    
        mcmcInfoIID = resample_chains_iid(mcmcInfoIID);                      
        % use Gibbs sampling to update hyperparameters    
        mcmcInfoIID = update_hmm_parameters_gibbs(mcmcInfoIID);
    end