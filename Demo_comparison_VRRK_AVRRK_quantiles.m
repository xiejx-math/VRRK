% This Matlab file compares VRRK-AS and AVRRK-AS on an inconsistent system.
% All problem, iteration, algorithm, and plotting parameters can be changed
% manually in the parameter section below.

clear;
close all;
clc;

% Locate this repository and make the shared implementations available.
repositoryRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(repositoryRoot,'src'));

%%% problem parameters
m=2000;
n=200;
rankA=200;
kappa=30;
delta=0.1;
solutionType="tail";        % Choose "random" or "tail".
seed=20260904;

%%% experiment parameters
maxIter=500*m;              % Total iterations used by each method.
runTime=20;                 % Number of independent random trials.
stopAtTolerance=false;      % Set true to stop when RSE <= tolerance.
tolerance=1e-12;
sampling="uniform";         % Choose "uniform" or "row-norm".

%%% automatic tuning parameters
autoTune=false;              % Set false to use the manual values below.
tuningSeed=seed+900000;     % Independent seed used only for parameter tuning.
tuningBudgets=[10 50 200];  % Matrix-access passes in successive halving.
tuningKeepFraction=1/3;
tuningStabilityPasses=500;

%%% manual VRRK-AS parameter, used only when autoTune=false
VRRK_cAlpha=1;

%%% manual AVRRK-AS parameters, used only when autoTune=false
AVRRK_cAlpha=4;
AVRRK_eta=0.25;
AVRRK_cTheta=1;
AVRRK_cGamma=1.25;
AVRRK_cBeta=1;
AVRRK_theta2=0.25;

%cAlpha=4, eta=0.25, cTheta=1, cGamma=1.25, cBeta=1, theta2=0.25

%%% plotting and output parameters
lineWidth=2.4;
fontSize=15;
saveFigures=true;
outputFolder=fullfile(repositoryRoot,'results','vrrk_avrrk');
outputFolder=new_result_folder(outputFolder,sprintf('%s_m%d_n%d_r%d_k%g_delta%g_seed%d',solutionType,m,n,rankA,kappa,delta,seed));
% Preserve editable settings and seeds before running any trials.
save(fullfile(outputFolder,'run_settings.mat'));

%%% check manually supplied parameters
if rankA>min(m,n)
    error('rankA must not exceed min(m,n).');
end
if maxIter<1 || maxIter~=round(maxIter)
    error('maxIter must be a positive integer.');
end
if ~ismember(solutionType,["random","tail"])
    error('solutionType must be "random" or "tail".');
end
if ~ismember(sampling,["uniform","row-norm"])
    error('sampling must be "uniform" or "row-norm".');
end
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end

%%% tune parameters or construct them from the manual values
if autoTune
    % Use the same matrix class and right-hand-side construction as the test.
    % The independent seed prevents tuning directly on a reported test trial.
    tuneProb=make_synthetic_problem(m,n,rankA,kappa,delta,...
        tuningSeed,solutionType);
    cfgTune=experiment_config("quick");
    cfgTune.seed=tuningSeed;
    cfgTune.sampling=sampling;
    cfgTune.tuning.budgets=tuningBudgets;
    cfgTune.tuning.keepFraction=tuningKeepFraction;
    cfgTune.tuning.stabilityPasses=tuningStabilityPasses;
    [parVRRK,parAVRRK,tuningLog]=tune_parameters(tuneProb,cfgTune);
    writetable(tuningLog,fullfile(outputFolder,'tuning_log.csv'));
else
    % Use the values specified in the manual parameter section above.
    parVRRK=struct('cAlpha',VRRK_cAlpha); %#ok<UNRCH>
    parAVRRK=struct('cAlpha',AVRRK_cAlpha,...
        'eta',AVRRK_eta,...
        'cTheta',AVRRK_cTheta,...
        'cGamma',AVRRK_cGamma,...
        'cBeta',AVRRK_cBeta,...
        'theta2',AVRRK_theta2);
    tuningLog=table();
end

%%% compute realized parameters on the tuning problem for reporting
if ~exist('tuneProb','var')
    tuneProb=make_synthetic_problem(m,n,rankA,kappa,delta,...
        tuningSeed,solutionType);
end
tuneP=select_row_probability(tuneProb,sampling);
tuneL=max(tuneProb.rowNorm2./max(tuneP,eps));
VRRK_alpha=parVRRK.cAlpha/tuneL;
AVRRK_alpha=parAVRRK.cAlpha/tuneL;
AVRRK_theta1=min(.25,parAVRRK.cTheta*sqrt(AVRRK_alpha));
AVRRK_theta2_real=min(parAVRRK.theta2,1-AVRRK_theta1);
AVRRK_gamma=parAVRRK.cGamma*AVRRK_alpha/2;
AVRRK_beta=min(1,max(0,1-parAVRRK.cBeta*AVRRK_gamma));

%%% print and save the selected dimensionless and realized parameters
fprintf('\nSelected parameters from the tuning instance\n');
fprintf('VRRK-AS:  cAlpha=%g, alpha=%.6e\n',...
    parVRRK.cAlpha,VRRK_alpha);
fprintf(['AVRRK-AS: cAlpha=%g, eta=%g, cTheta=%g, cGamma=%g, ',...
    'cBeta=%g, theta2=%g\n'],...
    parAVRRK.cAlpha,parAVRRK.eta,parAVRRK.cTheta,...
    parAVRRK.cGamma,parAVRRK.cBeta,parAVRRK.theta2);
fprintf('Realized: alpha=%.6e, theta1=%.6e, theta2=%.6e, gamma=%.6e, beta=%.6e\n',...
    AVRRK_alpha,AVRRK_theta1,AVRRK_theta2_real,...
    AVRRK_gamma,AVRRK_beta);

Method=["VRRK-AS";"AVRRK-AS"];
cAlpha=[parVRRK.cAlpha;parAVRRK.cAlpha];
eta=[NaN;parAVRRK.eta];
cTheta=[NaN;parAVRRK.cTheta];
cGamma=[NaN;parAVRRK.cGamma];
cBeta=[NaN;parAVRRK.cBeta];
theta2=[NaN;parAVRRK.theta2];
alpha=[VRRK_alpha;AVRRK_alpha];
theta1=[NaN;AVRRK_theta1];
gamma=[NaN;AVRRK_gamma];
beta=[NaN;AVRRK_beta];
parameterTable=table(Method,cAlpha,eta,cTheta,cGamma,cBeta,...
    theta2,alpha,theta1,gamma,beta);
writetable(parameterTable,fullfile(outputFolder,'selected_parameters.csv'));
save(fullfile(outputFolder,'selected_parameters.mat'),...
    'parVRRK','parAVRRK','parameterTable','tuningLog');

%%% allocate storage for histories and terminal statistics
H=struct('VRRK',{{}},'AVRRK',{{}});
VRRK_RSE=zeros(runTime,1);
AVRRK_RSE=zeros(runTime,1);
VRRK_Time=zeros(runTime,1);
AVRRK_Time=zeros(runTime,1);
VRRK_Iter=zeros(runTime,1);
AVRRK_Iter=zeros(runTime,1);

%%% run and store independent numerical trials
for jj=1:runTime
    % Generate a new matrix and inconsistent right-hand side for each trial.
    prob=make_synthetic_problem(m,n,rankA,kappa,delta,...
        seed+jj,solutionType);

    % Use separate budgets so that both methods perform maxIter iterations.
    % For this dense synthetic matrix, VRRK uses one row access per iteration
    % and AVRRK uses two row accesses per iteration.
    cfgVRRK=experiment_config("quick");
    cfgVRRK.sampling=sampling;
    cfgVRRK.maxPasses=maxIter/m;
    cfgVRRK.maxIterations=maxIter;
    cfgVRRK.recordPass=1;
    cfgVRRK.divergenceRSE=1e12;

    cfgAVRRK=cfgVRRK;
    cfgAVRRK.maxPasses=2*maxIter/m;
    cfgAVRRK.recordPass=2;

    if stopAtTolerance
        cfgVRRK.tol=tolerance; %#ok<UNRCH>
        cfgAVRRK.tol=tolerance;
    else
        % A zero tolerance forces both methods to use the full iteration budget.
        cfgVRRK.tol=0;
        cfgAVRRK.tol=0;
    end

    % Reset RNG before each method; both draw online from the same row stream.
    rng(seed+100000+jj,'twister');
    hVRRK=run_vrrk_as(prob,parVRRK,cfgVRRK);
    rng(seed+100000+jj,'twister');
    hAVRRK=run_avrrk_as(prob,parAVRRK,cfgAVRRK);

    % Use the same epoch definition as in the rank-deficient comparison.
    % Here k denotes the iteration counter of the corresponding method.
    hVRRK.epoch=hVRRK.iter/m;
    hAVRRK.epoch=hAVRRK.iter/m;

    % Store complete histories for median epoch and time curves.
    H.VRRK{end+1}=hVRRK;
    H.AVRRK{end+1}=hAVRRK;

    % Store terminal results for the printed summary.
    VRRK_RSE(jj)=hVRRK.rse(end);
    AVRRK_RSE(jj)=hAVRRK.rse(end);
    VRRK_Time(jj)=hVRRK.time(end);
    AVRRK_Time(jj)=hAVRRK.time(end);
    VRRK_Iter(jj)=hVRRK.iter(end);
    AVRRK_Iter(jj)=hAVRRK.iter(end);
end

% Save complete trajectories so every pointwise band can be verified later.
save(fullfile(outputFolder,'all_trial_histories.mat'),'H');

%%% print the median result over all trials
fprintf('\nVRRK-AS versus AVRRK-AS\n');
fprintf('m=%d, n=%d, rank=%d, kappa=%g, delta=%g, type=%s\n',...
    m,n,rankA,kappa,delta,solutionType);
fprintf('Method       Iterations       Time(s)      final RSE\n');
fprintf('VRRK-AS      %10.1f    %10.4f    %10.4e\n',...
    median(VRRK_Iter),median(VRRK_Time),median(VRRK_RSE));
fprintf('AVRRK-AS     %10.1f    %10.4f    %10.4e\n',...
    median(AVRRK_Iter),median(AVRRK_Time),median(AVRRK_RSE));

%%% define plotting data and colors in the same form as the other demos
methodNames={'VRRK-AS','AVRRK-AS'};
methodRuns={H.VRRK,H.AVRRK};
lineColors=[0.0000 0.4470 0.7410;
            0.8500 0.3250 0.0980];
lightColors=0.78+0.22*lineColors;
darkColors=0.48+0.52*lineColors;

%%% plot RSE versus epoch using the common multi-method plotting function
figEpoch=figure('Color','w');
axEpoch=axes(figEpoch);
hold(axEpoch,'on');
epochHandles=plot_multi_epoch_quantiles(axEpoch,methodRuns,...
    lineColors,lightColors,darkColors,lineWidth);
set(axEpoch,'YScale','log');
grid(axEpoch,'on');
box(axEpoch,'on');
xlabel(axEpoch,'Epochs $(k/m)$','Interpreter','latex');
ylabel(axEpoch,'RSE','Interpreter','latex');
title(axEpoch,sprintf('$%s\\ solution,\\ \\kappa=%g$',...
    solutionType,kappa),'Interpreter','latex');
legend(axEpoch,epochHandles,methodNames,...
    'Location','best','Interpreter','latex');
set(axEpoch,'FontSize',fontSize);
if saveFigures
    exportgraphics(figEpoch,fullfile(outputFolder,...
        sprintf('Demo_VRRK_AVRRK_quantiles_%s_k%g_epoch.png',...
        solutionType,kappa)),'Resolution',300);
end

%%% plot time-to-accuracy quantiles using the common plotting function
figTime=figure('Color','w');
axTime=axes(figTime);
hold(axTime,'on');
timeHandles=plot_multi_time_quantiles(axTime,methodRuns,...
    lineColors,lightColors,darkColors,lineWidth);
set(axTime,'YScale','log');
grid(axTime,'on');
box(axTime,'on');
xlabel(axTime,'Elapsed time (seconds)','Interpreter','latex');
ylabel(axTime,'RSE','Interpreter','latex');
title(axTime,sprintf('$%s\\ solution,\\ \\kappa=%g$',...
    solutionType,kappa),'Interpreter','latex');
legend(axTime,timeHandles,methodNames,...
    'Location','best','Interpreter','latex');
set(axTime,'FontSize',fontSize);
if saveFigures
    exportgraphics(figTime,fullfile(outputFolder,...
        sprintf('Demo_VRRK_AVRRK_quantiles_%s_k%g_time.png',...
        solutionType,kappa)),'Resolution',300);
end
function p = select_row_probability(prob,sampling)
%SELECT_ROW_PROBABILITY Return the manually requested row distribution.
if sampling=="row-norm"
    p=prob.rowP;
else
    p=prob.uniformP;
end
end
