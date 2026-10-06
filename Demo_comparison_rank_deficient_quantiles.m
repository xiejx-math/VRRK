% This Matlab file compares VRRK-AS, AVRRK-AS, and REK on synthetic
% rank-deficient inconsistent linear systems.

clear;
close all;
clc;

% Locate this repository and make the shared implementations available.
repositoryRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(repositoryRoot,'src'));

%%% problem parameters
m=2000;
n=200;
rankValues=[150 100 50];
kappa=100;
delta=0.1;
solutionType="random";           % Choose "random" or "tail".
seed=20260909;

%%% experiment parameters
maxEpoch=500;
rekBudgetMultiplier=8;           % Longer safeguard for the slower baseline.
runTime=20;
stopAtTolerance=true;
tolerance=1e-12;
sampling="uniform";              % Sampling for VRRK-AS and AVRRK-AS.

%%% automatic parameter tuning for VRRK-AS and AVRRK-AS
autoTune=true;
tuningBudgets=[10 50 200];
tuningKeepFraction=1/3;
tuningStabilityPasses=500;

%%% manual parameters used only when autoTune=false
VRRK_cAlpha=1;
AVRRK_cAlpha=2;
AVRRK_eta=0.5;
AVRRK_cTheta=1;
AVRRK_cGamma=1;
AVRRK_cBeta=1;
AVRRK_theta2=0.5;

%%% plotting and output parameters
lineWidth=2.4;
fontSize=15;
saveFigures=true;
outputFolder=fullfile(repositoryRoot,'results','rank_deficient');
outputFolder=new_result_folder(outputFolder,sprintf('%s_m%d_n%d_k%g_delta%g_seed%d',solutionType,m,n,kappa,delta,seed));
% Preserve editable settings and seeds before running any trials.
save(fullfile(outputFolder,'run_settings.mat'));

%%% validate manually supplied settings
if any(rankValues<1) || any(rankValues>=n)
    error('Every rank must satisfy 1 <= rank < n.');
end
if m<n
    error('This demo assumes an overdetermined system with m >= n.');
end
if kappa<1
    error('kappa must be at least one.');
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

%%% define consistent colors for all figures
methodNames={'VRRK-AS','AVRRK-AS','REK'};
lineColors=[0.0000 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.4660 0.6740 0.1880];
lightColors=0.78+0.22*lineColors;
darkColors=0.48+0.52*lineColors;

%%% allocate the summary table
summaryRows=struct('rankA',{},'kappa',{},'method',{},...
    'medianRSE',{},'q25RSE',{},'q75RSE',{},...
    'medianEpoch',{},'medianTime',{});
summaryIndex=0;

%%% run the experiment for every requested numerical rank
for rankIndex=1:numel(rankValues)
    rankA=rankValues(rankIndex);
    caseSeed=seed+1000000*rankIndex;
    tuningSeed=caseSeed+900000;
    caseFolder=fullfile(outputFolder,sprintf('rank_%d',rankA));
    if ~exist(caseFolder,'dir')
        mkdir(caseFolder);
    end

    fprintf('\n============================================================\n');
    fprintf(['Rank-deficient comparison: m=%d, n=%d, r=%d, ',...
        'kappa=%g, delta=%g, type=%s\n'],...
        m,n,rankA,kappa,delta,solutionType);

    % Generate an independent problem instance used only for tuning.
    tuneProb=make_synthetic_problem(m,n,rankA,kappa,delta,...
        tuningSeed,solutionType);

    % Tune the two proposed algorithms independently of reported trials.
    if autoTune
        cfgTune=experiment_config("quick");
        cfgTune.seed=tuningSeed;
        cfgTune.sampling=sampling;
        cfgTune.tuning.budgets=tuningBudgets;
        cfgTune.tuning.keepFraction=tuningKeepFraction;
        cfgTune.tuning.stabilityPasses=tuningStabilityPasses;
        [parVRRK,parAVRRK,tuningLog]=tune_parameters(tuneProb,cfgTune);
        writetable(tuningLog,fullfile(caseFolder,'tuning_log.csv'));
    else
        parVRRK=struct('cAlpha',VRRK_cAlpha); %#ok<UNRCH>
        parAVRRK=struct('cAlpha',AVRRK_cAlpha,...
            'eta',AVRRK_eta,...
            'cTheta',AVRRK_cTheta,...
            'cGamma',AVRRK_cGamma,...
            'cBeta',AVRRK_cBeta,...
            'theta2',AVRRK_theta2);
        tuningLog=table();
    end

    % Convert the dimensionless values into realized algorithm parameters.
    tuneP=select_rank_deficient_row_probability(tuneProb,sampling);
    tuneL=max(tuneProb.rowNorm2./max(tuneP,eps));
    VRRK_alpha=parVRRK.cAlpha/tuneL;
    AVRRK_alpha=parAVRRK.cAlpha/tuneL;
    AVRRK_theta1=min(0.25,parAVRRK.cTheta*sqrt(AVRRK_alpha));
    AVRRK_theta2_real=min(parAVRRK.theta2,1-AVRRK_theta1);
    AVRRK_gamma=parAVRRK.cGamma*AVRRK_alpha/2;
    AVRRK_beta=min(1,max(0,1-parAVRRK.cBeta*AVRRK_gamma));

    fprintf('VRRK-AS: cAlpha=%g, alpha=%.6e\n',...
        parVRRK.cAlpha,VRRK_alpha);
    fprintf(['AVRRK-AS: cAlpha=%g, eta=%g, cTheta=%g, ',...
        'cGamma=%g, cBeta=%g, theta2=%g\n'],...
        parAVRRK.cAlpha,parAVRRK.eta,parAVRRK.cTheta,...
        parAVRRK.cGamma,parAVRRK.cBeta,parAVRRK.theta2);
    fprintf(['Realized AVRRK-AS: alpha=%.6e, theta1=%.6e, ',...
        'theta2=%.6e, gamma=%.6e, beta=%.6e\n'],...
        AVRRK_alpha,AVRRK_theta1,AVRRK_theta2_real,...
        AVRRK_gamma,AVRRK_beta);

    % Save dimensionless and realized parameters for reproducibility.
    Method=["VRRK-AS";"AVRRK-AS"];
    cAlpha=[parVRRK.cAlpha;parAVRRK.cAlpha];
    eta=[NaN;parAVRRK.eta];
    cTheta=[NaN;parAVRRK.cTheta];
    cGamma=[NaN;parAVRRK.cGamma];
    cBeta=[NaN;parAVRRK.cBeta];
    theta2=[NaN;parAVRRK.theta2];
    alpha=[VRRK_alpha;AVRRK_alpha];
    theta1=[NaN;AVRRK_theta1];
    theta2Realized=[NaN;AVRRK_theta2_real];
    gamma=[NaN;AVRRK_gamma];
    beta=[NaN;AVRRK_beta];
    parameterTable=table(Method,cAlpha,eta,cTheta,cGamma,cBeta,theta2,...
        alpha,theta1,theta2Realized,gamma,beta);
    writetable(parameterTable,fullfile(caseFolder,'selected_parameters.csv'));

    % Allocate complete histories and terminal quantities.
    H=struct('VRRK',{{}},'AVRRK',{{}},'REK',{{}});
    finalRSE=zeros(runTime,3);
    finalEpoch=zeros(runTime,3);
    finalTime=zeros(runTime,3);

    for trial=1:runTime
        fprintf('  trial %2d of %2d\n',trial,runTime);

        % Generate an independent rank-deficient inconsistent system.
        prob=make_synthetic_problem(m,n,rankA,kappa,delta,...
            caseSeed+trial,solutionType);

        % VRRK-AS uses one row inner product in every iteration.
        cfgVRRK=experiment_config("quick");
        cfgVRRK.sampling=sampling;
        cfgVRRK.maxPasses=maxEpoch;
        cfgVRRK.recordPass=1;
        cfgVRRK.divergenceRSE=1e12;

        % AVRRK-AS uses two row inner products in every iteration.
        cfgAVRRK=cfgVRRK;
        cfgAVRRK.maxPasses=2*maxEpoch;
        cfgAVRRK.recordPass=2;

        % REK performs one row update and one column update per iteration.
        cfgREK=cfgVRRK;
        cfgREK.maxPasses=rekBudgetMultiplier*maxEpoch;
        cfgREK.recordPass=1;

        if stopAtTolerance
            cfgVRRK.tol=tolerance;
            cfgAVRRK.tol=tolerance;
            cfgREK.tol=tolerance;
        else
            cfgVRRK.tol=0; %#ok<UNRCH>
            cfgAVRRK.tol=0;
            cfgREK.tol=0;
        end

        % Couple the row-index sequences of VRRK-AS and AVRRK-AS.
        rng(caseSeed+100000+trial,'twister');
        hVRRK=run_vrrk_as(prob,parVRRK,cfgVRRK);
        rng(caseSeed+100000+trial,'twister');
        hAVRRK=run_avrrk_as(prob,parAVRRK,cfgAVRRK);

        % Use an independent reproducible stream for REK.
        rng(caseSeed+200000+trial,'twister');
        hREK=run_rek(prob,cfgREK);

        % Apply the method-dependent epoch definitions from the manuscript.
        hVRRK.epoch=hVRRK.iter/m;
        hAVRRK.epoch=hAVRRK.iter/m;
        hREK.epoch=hREK.iter/m+hREK.iter/n;

        H.VRRK{end+1}=hVRRK;
        H.AVRRK{end+1}=hAVRRK;
        H.REK{end+1}=hREK;

        trialHistories={hVRRK,hAVRRK,hREK};
        for methodIndex=1:3
            history=trialHistories{methodIndex};
            finalRSE(trial,methodIndex)=history.rse(end);
            finalEpoch(trial,methodIndex)=history.epoch(end);
            finalTime(trial,methodIndex)=history.time(end);
        end
    end

    % Preserve complete trajectories for later verification and replotting.
    save(fullfile(caseFolder,'all_trial_histories.mat'),...
        'H','parVRRK','parAVRRK','tuningLog','parameterTable',...
        'm','n','rankA','kappa','delta','solutionType',...
        'maxEpoch','runTime');

    % Add terminal medians and quartiles to the combined summary.
    for methodIndex=1:3
        summaryIndex=summaryIndex+1;
        summaryRows(summaryIndex).rankA=rankA;
        summaryRows(summaryIndex).kappa=kappa;
        summaryRows(summaryIndex).method=string(methodNames{methodIndex});
        summaryRows(summaryIndex).medianRSE=median(finalRSE(:,methodIndex));
        summaryRows(summaryIndex).q25RSE=...
            sample_percentile(finalRSE(:,methodIndex),25);
        summaryRows(summaryIndex).q75RSE=...
            sample_percentile(finalRSE(:,methodIndex),75);
        summaryRows(summaryIndex).medianEpoch=median(finalEpoch(:,methodIndex));
        summaryRows(summaryIndex).medianTime=median(finalTime(:,methodIndex));
    end

    % Arrange histories in the same order as the plot legend.
    methodRuns={H.VRRK,H.AVRRK,H.REK};

    % Plot RSE versus the prescribed method-dependent epochs.
    figEpoch=figure('Color','w');
    axEpoch=axes(figEpoch); %#ok<LAXES>
    hold(axEpoch,'on');
    epochHandles=plot_multi_epoch_quantiles(axEpoch,methodRuns,...
        lineColors,lightColors,darkColors,lineWidth);
    set(axEpoch,'YScale','log');
    grid(axEpoch,'on');
    box(axEpoch,'on');
    xlabel(axEpoch,'Epochs','Interpreter','latex');
    ylabel(axEpoch,'RSE','Interpreter','latex');
    title(axEpoch,sprintf('$r=%d,\\ \\kappa=%g$',rankA,kappa),...
        'Interpreter','latex');
    legend(axEpoch,epochHandles,methodNames,...
        'Location','best','Interpreter','latex');
    set(axEpoch,'FontSize',fontSize);
    if saveFigures
        exportgraphics(figEpoch,fullfile(caseFolder,sprintf(...
            'rank_%d_%s_kappa%g_epochs.png',...
            rankA,solutionType,kappa)),'Resolution',300);
    end

    % Plot elapsed-time quantiles at fixed RSE levels.
    figTime=figure('Color','w');
    axTime=axes(figTime); %#ok<LAXES>
    hold(axTime,'on');
    timeHandles=plot_multi_time_quantiles(axTime,methodRuns,...
        lineColors,lightColors,darkColors,lineWidth);
    set(axTime,'YScale','log');
    grid(axTime,'on');
    box(axTime,'on');
    xlabel(axTime,'Elapsed time (seconds)','Interpreter','latex');
    ylabel(axTime,'RSE','Interpreter','latex');
    title(axTime,sprintf('$r=%d,\\ \\kappa=%g$',rankA,kappa),...
        'Interpreter','latex');
    legend(axTime,timeHandles,methodNames,...
        'Location','best','Interpreter','latex');
    set(axTime,'FontSize',fontSize);
    if saveFigures
        exportgraphics(figTime,fullfile(caseFolder,sprintf(...
            'rank_%d_%s_kappa%g_time.png',...
            rankA,solutionType,kappa)),'Resolution',300);
    end
end

%%% write and display the combined numerical summary
summaryTable=struct2table(summaryRows);
writetable(summaryTable,fullfile(outputFolder,...
    'rank_deficient_summary.csv'));
disp(summaryTable);
fprintf('\nCompleted. Results are stored in %s.\n',outputFolder);

function p = select_rank_deficient_row_probability(prob,sampling)
%SELECT_RANK_DEFICIENT_ROW_PROBABILITY Return the requested distribution.
if sampling=="row-norm"
    p=prob.rowP;
else
    p=prob.uniformP;
end
end
