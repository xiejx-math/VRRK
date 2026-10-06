% This Matlab file compares VRRK-AS, AVRRK-AS, REK, and RCD on
% full-column-rank synthetic inconsistent linear systems.
% All principal experiment parameters can be changed in the section below.

clear;
close all;
clc;

% Locate this repository and make the shared implementations available.
repositoryRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(repositoryRoot,'src'));

%%% problem parameters
m=2000;
n=200;
rankA=n;                         % RCD is used only when rankA=n.
kappaValues=[ 300];
%kappaValues=[10 100 300];
delta=0.1;
solutionType="random";           % Choose "random" or "tail".
seed=20260907;

%%% experiment parameters
maxEpoch=1000;
baselineEpochMultiplier=5;       % Longer safeguards for REK and RCD.
runTime=20;
stopAtTolerance=false;
tolerance=1e-12;
sampling="uniform";              % Sampling for VRRK-AS and AVRRK-AS.

%%% automatic tuning parameters for VRRK-AS and AVRRK-AS
autoTune=true;
tuningBudgets=[10 50 200];
tuningKeepFraction=1/3;
tuningStabilityPasses=500;

%%% manual parameters used only when autoTune=false
VRRK_cAlpha=1;
AVRRK_cAlpha=4;
AVRRK_eta=0.25;
AVRRK_cTheta=1;
AVRRK_cGamma=1.25;
AVRRK_cBeta=1;
AVRRK_theta2=0.25;

% tail, kappa=30
%cAlpha=4, eta=0.25, cTheta=1, cGamma=1.25, cBeta=1, theta2=0.25

%%% plotting and output parameters
lineWidth=2.4;
fontSize=15;
saveFigures=true;
outputFolder=fullfile(repositoryRoot,'results','full_rank');
outputFolder=new_result_folder(outputFolder,sprintf('%s_m%d_n%d_delta%g_seed%d',solutionType,m,n,delta,seed));
% Preserve editable settings and seeds before running any trials.
save(fullfile(outputFolder,'run_settings.mat'));

%%% validate user-supplied settings
if rankA~=n || m<n
    error('This full-rank RCD experiment requires m>=n and rankA=n.');
end
if any(kappaValues<1)
    error('Every condition number must be at least one.');
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
methodNames={'VRRK-AS','AVRRK-AS','REK','RCD'};
lineColors=[0.0000 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.4660 0.6740 0.1880;
            0.4940 0.1840 0.5560];
lightColors=0.78+0.22*lineColors;
darkColors=0.48+0.52*lineColors;

%%% allocate the summary table
summaryRows=struct('kappa',{},'method',{},'medianRSE',{},...
    'q25RSE',{},'q75RSE',{},'medianEpoch',{},'medianTime',{});
summaryIndex=0;

%%% run the complete experiment for every requested condition number
for kk=1:numel(kappaValues)
    kappa=kappaValues(kk);
    caseSeed=seed+1000000*kk;
    tuningSeed=caseSeed+900000;
    caseFolder=fullfile(outputFolder,sprintf('kappa_%g',kappa));
    if ~exist(caseFolder,'dir')
        mkdir(caseFolder);
    end

    fprintf('\n============================================================\n');
    fprintf('Full-rank comparison: kappa=%g, delta=%g, type=%s\n',...
        kappa,delta,solutionType);

    % Tune only the proposed methods on an independent problem instance.
    if autoTune
        tuneProb=make_synthetic_problem(m,n,rankA,kappa,delta,...
            tuningSeed,solutionType);
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
        tuneProb=make_synthetic_problem(m,n,rankA,kappa,delta,...
            tuningSeed,solutionType);
    end

    fprintf('VRRK-AS: cAlpha=%g\n',parVRRK.cAlpha);
    fprintf(['AVRRK-AS: cAlpha=%g, eta=%g, cTheta=%g, ',...
        'cGamma=%g, cBeta=%g, theta2=%g\n'],...
        parAVRRK.cAlpha,parAVRRK.eta,parAVRRK.cTheta,...
        parAVRRK.cGamma,parAVRRK.cBeta,parAVRRK.theta2);

    % Convert the selected dimensionless values on the tuning instance.
    tuneP=select_demo_row_probability(tuneProb,sampling);
    tuneL=max(tuneProb.rowNorm2./max(tuneP,eps));
    VRRK_alpha=parVRRK.cAlpha/tuneL;
    AVRRK_alpha=parAVRRK.cAlpha/tuneL;
    AVRRK_theta1=min(0.25,parAVRRK.cTheta*sqrt(AVRRK_alpha));
    AVRRK_theta2_real=min(parAVRRK.theta2,1-AVRRK_theta1);
    AVRRK_gamma=parAVRRK.cGamma*AVRRK_alpha/2;
    AVRRK_beta=min(1,max(0,1-parAVRRK.cBeta*AVRRK_gamma));
    fprintf('VRRK-AS realized alpha=%.6e\n',VRRK_alpha);
    fprintf(['AVRRK-AS realized alpha=%.6e, theta1=%.6e, ',...
        'theta2=%.6e, gamma=%.6e, beta=%.6e\n'],...
        AVRRK_alpha,AVRRK_theta1,AVRRK_theta2_real,...
        AVRRK_gamma,AVRRK_beta);

    % Save both dimensionless and realized parameters for reproducibility.
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

    % Allocate complete trajectory histories for all four methods.
    H=struct('VRRK',{{}},'AVRRK',{{}},'REK',{{}},'RCD',{{}});

    % Store terminal quantities used in the numerical summary.
    finalRSE=zeros(runTime,4);
    finalEpoch=zeros(runTime,4);
    finalTime=zeros(runTime,4);

    for trial=1:runTime
        fprintf('  trial %2d of %2d\n',trial,runTime);

        % Every trial uses an independently generated A and b.
        prob=make_synthetic_problem(m,n,rankA,kappa,delta,...
            caseSeed+trial,solutionType);

        % A VRRK-AS iteration costs one row access and has epoch k/m.
        cfgVRRK=experiment_config("quick");
        cfgVRRK.sampling=sampling;
        cfgVRRK.maxPasses=maxEpoch;
        cfgVRRK.recordPass=1;
        cfgVRRK.divergenceRSE=1e12;

        % AVRRK-AS uses two row inner products per iteration. Doubling the
        % access budget allows it to perform the same maximum k/m epochs.
        cfgAVRRK=cfgVRRK;
        cfgAVRRK.maxPasses=2*maxEpoch;
        cfgAVRRK.recordPass=2;

        % Give the slower baselines a longer terminal safeguard. Their
        % displayed horizontal coordinates still use the epoch definitions.
        cfgRCD=cfgVRRK;
        cfgRCD.maxPasses=baselineEpochMultiplier*maxEpoch;
        cfgRCD.recordPass=1;

        cfgREK=cfgVRRK;
        cfgREK.maxPasses=baselineEpochMultiplier*maxEpoch;
        cfgREK.recordPass=1;

        if stopAtTolerance
            cfgVRRK.tol=tolerance; %#ok<UNRCH>
            cfgAVRRK.tol=tolerance;
            cfgREK.tol=tolerance;
            cfgRCD.tol=tolerance;
        else
            cfgVRRK.tol=0;
            cfgAVRRK.tol=0;
            cfgREK.tol=0;
            cfgRCD.tol=0;
        end

        % Use the same row-index sequence for the two proposed methods.
        rng(caseSeed+100000+trial,'twister');
        hVRRK=run_vrrk_as(prob,parVRRK,cfgVRRK);
        rng(caseSeed+100000+trial,'twister');
        hAVRRK=run_avrrk_as(prob,parAVRRK,cfgAVRRK);

        % Use reproducible independent streams for REK and RCD.
        rng(caseSeed+200000+trial,'twister');
        hREK=run_rek(prob,cfgREK);
        rng(caseSeed+300000+trial,'twister');
        hRCD=run_rcd(prob,cfgRCD);

        % Replace the generic row epoch by the method-specific definitions.
        hVRRK.epoch=hVRRK.iter/m;
        hAVRRK.epoch=hAVRRK.iter/m;
        hREK.epoch=hREK.iter/m+hREK.iter/n;
        hRCD.epoch=hRCD.iter/n;

        H.VRRK{end+1}=hVRRK;
        H.AVRRK{end+1}=hAVRRK;
        H.REK{end+1}=hREK;
        H.RCD{end+1}=hRCD;

        trialHistories={hVRRK,hAVRRK,hREK,hRCD};
        for methodIndex=1:4
            history=trialHistories{methodIndex};
            finalRSE(trial,methodIndex)=history.rse(end);
            finalEpoch(trial,methodIndex)=history.epoch(end);
            finalTime(trial,methodIndex)=history.time(end);
        end
    end

    save(fullfile(caseFolder,'all_trial_histories.mat'),...
        'H','parVRRK','parAVRRK','tuningLog','kappa','delta',...
        'solutionType','maxEpoch','runTime');

    % Append terminal medians and quartiles to the global summary.
    for methodIndex=1:4
        summaryIndex=summaryIndex+1;
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

    methodRuns={H.VRRK,H.AVRRK,H.REK,H.RCD};

    % Plot RSE statistics at common method-dependent epoch values.
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
    title(axEpoch,sprintf('$%s\\ solution,\\ \\kappa=%g$',...
        solutionType,kappa),'Interpreter','latex');
    legend(axEpoch,epochHandles,methodNames,...
        'Location','best','Interpreter','latex');
    set(axEpoch,'FontSize',fontSize);
    if saveFigures
        exportgraphics(figEpoch,fullfile(caseFolder,sprintf(...
            'full_rank_%s_kappa%g_epochs.png',solutionType,kappa)),...
            'Resolution',300);
    end

    % Plot median elapsed time versus median RSE at shared epoch checkpoints.
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
    title(axTime,sprintf('$%s\\ solution,\\ \\kappa=%g$',...
        solutionType,kappa),'Interpreter','latex');
    legend(axTime,timeHandles,methodNames,...
        'Location','best','Interpreter','latex');
    set(axTime,'FontSize',fontSize);
    if saveFigures
        exportgraphics(figTime,fullfile(caseFolder,sprintf(...
            'full_rank_%s_kappa%g_time.png',solutionType,kappa)),...
            'Resolution',300);
    end
end

%%% write the combined numerical summary
summaryTable=struct2table(summaryRows);
writetable(summaryTable,fullfile(outputFolder,'full_rank_summary.csv'));
disp(summaryTable);
fprintf('\nCompleted. Results are stored in %s.\n',outputFolder);

function p = select_demo_row_probability(prob,sampling)
%SELECT_DEMO_ROW_PROBABILITY Return the requested row distribution.
if sampling=="row-norm"
    p=prob.rowP;
else
    p=prob.uniformP;
end
end
