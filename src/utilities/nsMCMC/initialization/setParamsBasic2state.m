function trueParams = setParamsBasic2state

trueParams = struct;

trueParams.tres = 10;

% define parameters to be inferred
trueParams.R = [-.02, .05 ; .02 -.05]/2;
trueParams.A = expm(trueParams.R*trueParams.tres);
trueParams.v = [0, 4]' * trueParams.tres / 20;
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps = 8.3;
% mcmcInfo.trueParams = trueParams;

% add known hyperparameters
trueParams.nStates = size(trueParams.A,1);
trueParams.alpha_frac = 1.05/6;


