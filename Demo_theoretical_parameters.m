% Small validation run using Theorem 4.1 parameters (no parameter tuning).
clear;
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
prob=make_synthetic_problem(80,12,12,10,.1,20260913,"random");
cfg=experiment_config("quick");
cfg.maxPasses=20;
cfg.recordPass=1;
cfg.tol=0;
parameters=struct('mode',"theory");
values=avrrk_parameters(prob,parameters,cfg.sampling);
disp(values);
rng(20260914,'twister');
history=run_avrrk_as(prob,parameters,cfg);
fprintf('Theoretical parameters: final RSE %.6e after %d iterations.\n',...
    history.rse(end),history.iter(end));
% The conservative theoretical stepsize may lead to slow practical convergence.
