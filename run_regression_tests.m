function run_regression_tests()
%RUN_REGRESSION_TESTS Check scientific semantics, not just finite trajectories.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
% Preserve aligned checkpoints and never extrapolate missing leading data.
sample=struct('epoch',[2 3 5],'rse',[1 .1 .01],'time',[0 1 2]);
[gridCheck,~]=curve_statistics({sample,sample},'epoch');
assert(isequal(gridCheck,[2 3 5]));
later=struct('epoch',[3 4 6],'rse',[1 .2 .02],'time',[0 1 2]);
[gridCheck,~]=curve_statistics({sample,later},'epoch');
assert(gridCheck(1)==3 && gridCheck(end)==5);
caught=false;
try
    make_synthetic_problem(4,4,4,5,.1,17,"random");
catch exception
    caught=strcmp(exception.identifier,'VRRK:NoResidualSpace');
end
assert(caught,'A full-row-rank inconsistent construction must be rejected.');

% Repeated levels must select the first recorded attainment, never the last.
h=struct('time',[0 1 2 3],'rse',[1 .1 .1 .01],'epoch',[0 1 2 3]);
[~,stats]=time_curve_statistics({h},[1 .1 .01 .001]);
assert(isequal(stats.timeAtRSE,[0 1 3 Inf]));
h.rse=ones(1,4);
[~,stats]=time_curve_statistics({h});
assert(all(stats.median==0));
h=struct('time',0,'rse',1,'epoch',0);
[~,stats]=time_curve_statistics({h});
assert(all(stats.median==0));

% With 20 trials, the arithmetic median differs from the geometric median.
% Time curves keep the median below a poorly converging trial's final RSE.
timeRuns=cell(20,1);
for j=1:20
    value=.008;
    if j==20
        value=.03;
    end
    timeRuns{j}=struct('time',[0 1 2],'rse',[1 .1 value],'epoch',[0 1 2]);
end
[t,s]=curve_statistics(timeRuns,'time');
assert(t(end)==2 && abs(s.median(end)-.008)<1e-14);
assert(s.minimum(end)==.008 && s.maximum(end)==.03);
assert(s.q25(end)==.008 && s.q75(end)==.008);
timeRuns{20}.time=[0 .5 1];
[t,s]=curve_statistics(timeRuns,'time');
assert(t(end)==1 && abs(s.median(end)-.1)<1e-14);
% Plot objects must use time horizontally and RSE vertically, on a log axis.
[t,s,e]=median_time_rse_statistics(timeRuns);
[epochGrid,epochStats]=curve_statistics(timeRuns,'epoch');
assert(isequal(e,epochGrid) && isequal(s,epochStats));
assert(t(end)==2 && s.median(end)==.008 && all(diff(t)>=0));
shortRuns=timeRuns;
shortRuns{20}=struct('epoch',[0 1],'time',[0 .5],'rse',[1 .03]);
[~,shortStats,shortEpoch]=median_time_rse_statistics(shortRuns);
[expectedEpoch,expectedStats]=curve_statistics(shortRuns,'epoch');
assert(shortEpoch(end)==1 && isequal(shortEpoch,expectedEpoch));
assert(isequal(shortStats,expectedStats));
fig=figure('Visible','off');
cleanup=onCleanup(@() close(fig));
ax=axes(fig);
handles=plot_multi_time_quantiles(ax,{timeRuns},[0 0 1],[.8 .8 1],[.4 .4 1],2);
assert(strcmp(ax.YScale,'log'));
assert(isequal(handles.XData,t) && isequal(handles.YData,s.median));
assert(numel(findall(ax,'Type','patch'))==2);
clear cleanup;
runs=cell(20,1);
for j=1:20
    value=1e-4;
    if j>10
        value=1e-2;
    end
    runs{j}=struct('epoch',[0 1],'rse',[value value]);
end
[~,stats]=curve_statistics(runs,'epoch');
assert(max(abs(stats.median-0.00505))<1e-14);
assert(abs(sample_percentile(1:20,25)-5.5)<1e-14);
assert(abs(sample_percentile(1:20,75)-15.5)<1e-14);
if exist('prctile','file')
    assert(abs(sample_percentile(1:20,25)-prctile(1:20,25))<1e-14);
end

% The loader must transpose off-diagonal coordinates when completing symmetry.
A=mmread_simple(fullfile(root,'tests','symmetric_fixture.mtx'));
assert(isequal(full(A),[2 3 0;3 4 5;0 5 0]));

% A terminal state must be saved even before the first recording checkpoint.
prob=make_synthetic_problem(40,8,8,5,.1,17,"random");
cfg=experiment_config("quick");
cfg.tol=0;
cfg.maxPasses=100;
cfg.recordPass=100;
cfg.maxIterations=7;
cfg.divergenceRSE=1e12;
v=struct('cAlpha',.1);
a=struct('cAlpha',.1,'eta',.5,'cTheta',1,'cGamma',1,...
    'cBeta',1,'theta2',.5);
rng(19);
histories={run_vrrk_as(prob,v,cfg),run_avrrk_as(prob,a,cfg),...
    run_rek(prob,cfg),run_rcd(prob,cfg)};
for j=1:4
    result=histories{j};
    assert(result.iter(end)==7 && numel(result.iter)==2);
    expected=norm(result.xFinal-prob.xStar)^2/norm(prob.xStar)^2;
    assert(abs(result.rse(end)-expected)<1e-12);
    assert(result.stopReason=="budget" && ~result.failed);
    assert(result.samplingIncluded && all(diff(result.time)>=0));
end
assert(histories{3}.epoch(end)==7/prob.m+7/prob.n);
assert(histories{4}.epoch(end)==7/prob.n);

% Stopping is checked at diagnostics, and tol=0 disables convergence stopping.
cfg.recordPass=.01;
cfg.maxIterations=100;
cfg.tol=2;
result=run_vrrk_as(prob,v,cfg);
assert(result.stopReason=="tolerance" && result.iter(end)<100);
cfg.tol=0;
cfg.divergenceRSE=0;
result=run_rek(prob,cfg);
assert(result.failed && result.stopReason=="diverged");

% Theoretical parameters must use the measured smallest positive singular value.
theory=avrrk_parameters(prob,struct('mode',"theory"),"uniform");
assert(abs(theory.mu-1/25)<1e-12);
assert(abs(theory.beta-(1-theory.gamma*theory.mu))<1e-14);
assert(abs(theory.gamma-theory.alpha/(2*(theory.theta1+theory.alpha*theory.mu)))<1e-14);

% Stability failure must not silently return an unvalidated parameter set.
cfg=experiment_config("quick");
cfg.tuning.vrrkCAlpha=.1;
cfg.tuning.avStepScale=1;
cfg.tuning.eta=.5;
cfg.tuning.cTheta=1;
cfg.tuning.cGamma=1;
cfg.tuning.theta2=.5;
cfg.tuning.budgets=.1;
cfg.tuning.stabilityPasses=.2;
cfg.recordPass=.1;
cfg.divergenceRSE=0;
caught=false;
try
    tune_parameters(prob,cfg);
catch exception
    caught=strcmp(exception.identifier,'VRRK:TuningFailed');
end
assert(caught,'Tuning failure did not propagate.');

% Successful tuning includes both methods' stability evaluations and parameters.
cfg.divergenceRSE=1e12;
[~,~,log]=tune_parameters(prob,cfg);
assert(all(ismember({'cAlpha','seed','stage','failed'},log.Properties.VariableNames)));
assert(any(log.stage=="stability" & log.method=="VRRK"));
assert(any(log.stage=="stability" & log.method=="AVRRK"));
fprintf('All regression tests passed.\n');
% An excessive step must trigger bounded recovery, not return a bad run.
cfg.tuning.vrrkCAlpha=64;
cfg.tuning.budgets=[1 10];
cfg.tuning.stabilityPasses=20;
cfg.maxIterations=1;
[bestV,bestA,log]=tune_parameters(prob,cfg);
assert(bestV.cAlpha<64);
assert(any(contains(log.stage,"_recovery")));
assert(max(log.iterations)>1,'Demo iteration cap leaked into tuning.');
for method=["VRRK" "AVRRK"]
    chosen=bestV;
    if method=="AVRRK"
        chosen=bestA;
    end
    rows=log.stage=="stability" & log.method==method & log.cAlpha==chosen.cAlpha;
    assert(nnz(rows)==2 && all(log.finalRSE(rows)<1) && ~any(log.failed(rows)));
end
fprintf('Tuning recovery tests passed.\n');
% Exercise AVRRK recovery separately with an aggressive acceleration step.
cfg.tuning.vrrkCAlpha=.1;
cfg.tuning.avStepScale=640;
[~,~,log]=tune_parameters(prob,cfg);
assert(any(contains(log.stage,"_recovery") & log.method=="AVRRK"));
fprintf('AVRRK recovery test passed.\n');
end
