function trueParams = setParamsBasic3state(tres)

% mcmcInfo = struct;
trueParams = struct;
trueParams.discrete_data_flag = true;
trueParams.tres = tres;
trueParams.nStates = 3;

% define parameters to be inferred
trueParams.R = [-0.0200    0.0210         0;
                 0.0200   -0.0330    0.0700;
                      0    0.0120   -0.0700]; %NL: from appendix table S2 in Lammers 2020
                    
trueParams.A = expm(trueParams.R*trueParams.tres);
trueParams.v = [0, 2, 5]';
trueParams.sigma = 1.5;
[V, D] = eig(trueParams.A);
[~, mi] = max(real(diag(D)));
trueParams.pi0 = V(:,mi)/sum(V(:,mi));
trueParams.nSteps = 7;
trueParams.alpha_frac = 30/140;
