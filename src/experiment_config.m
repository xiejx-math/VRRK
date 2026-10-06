function cfg = experiment_config(mode)
%EXPERIMENT_CONFIG Create reproducible settings for all experiments.
if nargin==0
    mode="quick";
end

% Shared numerical, sampling, and output settings.
cfg.mode=string(mode);
cfg.seed=20260902;
cfg.delta=0.1;
cfg.tol=1e-8;
cfg.recordPass=1;
cfg.divergenceRSE=1e8;
cfg.sampling="uniform";
cfg.outputDir="results";

% Small defaults allow a fast end-to-end code check.
cfg.numTrials=3;
cfg.maxPasses=100;
cfg.synthetic.m=400;
cfg.synthetic.n=60;
cfg.synthetic.kappas=[1 10 30 100 300 1000];
cfg.synthetic.ranks=[60 45 30];
cfg.synthetic.deltas=[1e-2 1e-1 1];

% Successive-halving budgets and parameter grids.
cfg.tuning.budgets=[3 10 30];
cfg.tuning.keepFraction=1/3;
cfg.tuning.vrrkCAlpha=2.^(-10:2:0);
cfg.tuning.avStepScale=[1/4 1/2 1 2 4];
cfg.tuning.eta=[.25 .5 .75 .9 1];
cfg.tuning.cTheta=[.25 .5 1 2];
cfg.tuning.cGamma=[.5 .75 1 1.25 1.5];
cfg.tuning.theta2=[0 .25 .5 .75];
cfg.tuning.stabilityPasses=100;

% Fallback parameter values for direct method calls.
cfg.defaultVRRK=struct('cAlpha',1);
cfg.defaultAVRRK=struct('cAlpha',4,'eta',.5,'cTheta',.5,...
    'cGamma',1,'cBeta',1,'theta2',.5);

% Replace quick defaults by publication-scale settings.
if mode=="paper"
    cfg.numTrials=20;
    cfg.maxPasses=500;
    cfg.synthetic.m=2000;
    cfg.synthetic.n=200;
    cfg.synthetic.kappas=[1 3 10 30 100 300 1000];
    cfg.synthetic.ranks=[200 150 100 50];
    cfg.synthetic.deltas=[1e-3 1e-2 1e-1 1 10];
    cfg.tuning.budgets=[10 50 200];
    cfg.tuning.stabilityPasses=500;
end
end
