function trueParams = setParamsBasic3state(upsample_factor)

% mcmcInfo = struct;
trueParams = struct;

trueParams.tres_data = 20;
trueParams.tres_true = trueParams.tres_data/upsample_factor;
trueParams.upsample_factor = upsample_factor;

trueParams.nStates = 3;
% define parameters to be inferred
trueParams.R = [-0.0200    0.0210         0;
                 0.0200   -0.0330    0.0700;
                      0    0.0120   -0.0700]; %NL: from appendix table S2 in Lammers 2020
trueParams.A = expm(trueParams.R*trueParams.tres_true);
trueParams.v = [0, 2, 4]' / upsample_Factor;
trueParams.sigma = 0.8;
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps_data = 8.05;
trueParams.nsteps_true = trueParams.nSteps_data * upsample_factor;
trueParams.alpha_frac = 1.05/6;
