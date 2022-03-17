function trueParams = setParamsBasic3state

% mcmcInfo = struct;
trueParams = struct;

trueParams.tres = 20;
trueParams.nStates = 3;
% define parameters to be inferred
trueParams.R = [-.02, .04, 0; .02 -.05 .08; 0 .01 -.08];
trueParams.A = expm(trueParams.R*trueParams.tres);
trueParams.v = [0, 2, 4]';
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps = 8.05;
trueParams.alpha_frac = 1.05/6;
