function mcmcInfo = setParamsBasic3state

mcmcInfo = struct;
trueParams = struct;

mcmcInfo.tres = 20;

% define parameters to be inferred
trueParams.R = [-.02, .04, 0; .02 -.05 .08; 0 .01 -.08];
trueParams.A = expm(trueParams.R*mcmcInfo.tres);
trueParams.v = [0, 2, 4]';
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps = 6.9;
mcmcInfo.trueParams = trueParams;

% add known hyperparameters
mcmcInfo.nStates = size(trueParams.A,1);
mcmcInfo.alpha_frac = 1.05/6;
% mcmcInfo.trueParams.alpha = mcmcInfo.alpha_frac*mcmcInfo.trueParams.nSteps;
mcmcInfo.eps = 1e-2; % NL: note that this is not currently used
mcmcInfo.n_reps = 1;

% testing flags
mcmcInfo.ensembleInferenceFlag = 0;
mcmcInfo.testResampling = 0;
mcmcInfo.MHResampling = 0;
mcmcInfo.update_increment = 1; % sets how often parameter values are recorded in inference arrays

% calculate MS2 convolution kernel
mcmcInfo.trueParams.alpha = mcmcInfo.alpha_frac*mcmcInfo.trueParams.nSteps;
coeff_MS2 = ms2_loading_coeff_frac(mcmcInfo.trueParams.alpha, mcmcInfo.trueParams.nSteps, 12)';

% calculate scaled noise
f_mean = (mcmcInfo.trueParams.pi0'*mcmcInfo.trueParams.v)*sum(coeff_MS2);   
mcmcInfo.trueParams.sigma = f_mean/10;  