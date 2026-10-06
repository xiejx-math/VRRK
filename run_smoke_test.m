function run_smoke_test()
%RUN_SMOKE_TEST Verify the four solver implementations on a small problem.

repositoryRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(repositoryRoot,'src'));

seed=20260910;
prob=make_synthetic_problem(120,30,30,10,0.1,seed,"random");
cfg=experiment_config("quick");
cfg.maxPasses=3;
cfg.recordPass=1;
cfg.tol=0;
cfg.sampling="uniform";

% Couple the row samples used by the two proposed methods.
rng(seed+1,'twister');
indices=sample_discrete(prob.uniformP,ceil(cfg.maxPasses*prob.m));

hVRRK=run_vrrk_as(prob,cfg.defaultVRRK,cfg,indices);
hAVRRK=run_avrrk_as(prob,cfg.defaultAVRRK,cfg,indices);

% Use independent reproducible streams for the two baselines.
rng(seed+2,'twister');
hREK=run_rek(prob,cfg);
rng(seed+3,'twister');
hRCD=run_rcd(prob,cfg);

histories={hVRRK,hAVRRK,hREK,hRCD};
names={'VRRK-AS','AVRRK-AS','REK','RCD'};

fprintf('\nSmoke-test terminal RSE values\n');
for methodIndex=1:numel(histories)
    history=histories{methodIndex};
    assert(~isempty(history.rse) && all(isfinite(history.rse)),...
        '%s returned a nonfinite trajectory.',names{methodIndex});
    fprintf('  %-8s  %.4e\n',names{methodIndex},history.rse(end));
end

fprintf('Smoke test passed.\n');
end
