function trueParams = setParamsBasic2state(upsample_factor)

trueParams = struct;

trueParams.tres_data = 20;
trueParams.tres_true = trueParams.tres_data/upsample_factor;
trueParams.upsample_factor = upsample_factor;

% define parameters to be inferred
trueParams.R = [-.02, .05 ; .02 -.05];
trueParams.A = expm(trueParams.R*trueParams.tres_true);
trueParams.v = [0, 4]' / upsample_factor;
trueParams.sigma = 0.8;
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps_data = 8.05;
trueParams.nsteps_true = trueParams.nSteps_data * upsample_factor;

% add known hyperparameters
trueParams.nStates = size(trueParams.A,1);
trueParams.alpha_frac = 1.05/6;


