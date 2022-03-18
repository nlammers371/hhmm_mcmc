function trueParams = setParamsBasic3state

% mcmcInfo = struct;
trueParams = struct;

trueParams.tres = 20;
trueParams.nStates = 3;
% define parameters to be inferred
trueParams.R = [-0.0200    0.0210         0;
                 0.0200   -0.0330    0.0700;
                      0    0.0120   -0.0700]; %NL: from appendix table S2 in Lammers 2020
trueParams.A = expm(trueParams.R*trueParams.tres);
trueParams.v = [0, 2, 4]';
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps = 8.05;
trueParams.alpha_frac = 1.05/6;
