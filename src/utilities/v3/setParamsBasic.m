function mcmcInfo = setParamsBasic

mcmcInfo = struct;
trueParams = struct;

mcmcInfo.tres = 20;
% define parameters to be inferred
trueParams.R = [-.02, .04, 0; .02 -.05 .08; 0 .01 -.08];
trueParams.A = expm(trueParams.R*mcmcInfo.tres);
trueParams.v = [.05, 2, 4]';
trueParams.sigma = 0.5;
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));

mcmcInfo.trueParams = trueParams;

% add known hyperparameters
mcmcInfo.nStates = size(trueParams.A,1);
mcmcInfo.alpha = 1.05;
mcmcInfo.eps = 1e-2; % NL: note that this is not currently used
mcmcInfo.n_reps = 2;

mcmcInfo.testResampling = 0;
mcmcInfo.MHResampling = 0;