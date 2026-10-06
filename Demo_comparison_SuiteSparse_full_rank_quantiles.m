% This Matlab file compares VRRK-AS, AVRRK-AS, REK, and RCD on three
% full-column-rank matrices from the SuiteSparse Matrix Collection.

clear;
close all;
clc;

% Locate this repository and make the shared implementations available.
repositoryRoot=fileparts(mfilename('fullpath'));
addpath(fullfile(repositoryRoot,'src'));

%%% SuiteSparse matrices and inconsistent-system parameters
matrixNames=["ash958","well1033","illc1033"];
delta=0.1;
seed=20260909;
dataFolder=fullfile(repositoryRoot,'data','suitesparse');



%%% experiment parameters
maxEpoch=500;
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
outputFolder=fullfile(repositoryRoot,'results','suitesparse');
outputFolder=new_result_folder(outputFolder,sprintf('delta%g_seed%d',delta,seed));
% Preserve editable settings and seeds before running any trials.
save(fullfile(outputFolder,'run_settings.mat'));

%%% validate settings and locate the local matrix files
if ~ismember(sampling,["uniform","row-norm"])
    error('sampling must be "uniform" or "row-norm".');
end
if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end
matrixFiles=download_suitesparse(dataFolder,matrixNames);

%%% use consistent colors in all SuiteSparse figures
methodNames={'VRRK-AS','AVRRK-AS','REK','RCD'};
lineColors=[0.0000 0.4470 0.7410;
            0.8500 0.3250 0.0980;
            0.4660 0.6740 0.1880;
            0.4940 0.1840 0.5560];
lightColors=0.78+0.22*lineColors;
darkColors=0.48+0.52*lineColors;

%%% allocate the combined terminal-result table
summaryRows=struct('matrix',{},'m',{},'n',{},'nnzA',{},...
    'method',{},'medianRSE',{},'q25RSE',{},'q75RSE',{},...
    'medianEpoch',{},'medianTime',{});
summaryIndex=0;

%%% run all three full-column-rank SuiteSparse tests
for matrixIndex=1:numel(matrixNames)
    matrixName=matrixNames(matrixIndex);
    matrixFile=matrixFiles{matrixIndex};
    caseSeed=seed+1000000*matrixIndex;
    tuningSeed=caseSeed+900000;
    caseFolder=fullfile(outputFolder,matrixName);
    if ~exist(caseFolder,'dir')
        mkdir(caseFolder);
    end

    fprintf('\n============================================================\n');
    fprintf('SuiteSparse comparison: %s\n',matrixName);

    % Generate a separate right-hand side used only for parameter tuning.
    tuneProb=make_suitesparse_problem(matrixFile,delta,tuningSeed,true);

    % Tune the two proposed algorithms independently of all reported runs.
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

    % Convert the dimensionless parameters on the tuning problem.
    tuneP=select_suite_row_probability(tuneProb,sampling);
    tuneL=max(tuneProb.rowNorm2./max(tuneP,eps));
    VRRK_alpha=parVRRK.cAlpha/tuneL;
    AVRRK_alpha=parAVRRK.cAlpha/tuneL;
    AVRRK_theta1=min(0.25,parAVRRK.cTheta*sqrt(AVRRK_alpha));
    AVRRK_theta2_real=min(parAVRRK.theta2,1-AVRRK_theta1);
    AVRRK_gamma=parAVRRK.cGamma*AVRRK_alpha/2;
    AVRRK_beta=min(1,max(0,1-parAVRRK.cBeta*AVRRK_gamma));

    fprintf('Matrix size: m=%d, n=%d, nnz=%d\n',...
        tuneProb.m,tuneProb.n,tuneProb.nnz);
    fprintf('VRRK-AS: cAlpha=%g, alpha=%.6e\n',...
        parVRRK.cAlpha,VRRK_alpha);
    fprintf(['AVRRK-AS: cAlpha=%g, eta=%g, cTheta=%g, ',...
        'cGamma=%g, cBeta=%g, theta2=%g\n'],...
        parAVRRK.cAlpha,parAVRRK.eta,parAVRRK.cTheta,...
        parAVRRK.cGamma,parAVRRK.cBeta,parAVRRK.theta2);

    % Save selected and realized parameters for this matrix.
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

    % Allocate complete histories and terminal values for 20 trials.
    H=struct('VRRK',{{}},'AVRRK',{{}},'REK',{{}},'RCD',{{}});
    finalRSE=zeros(runTime,4);
    finalEpoch=zeros(runTime,4);
    finalTime=zeros(runTime,4);

    for trial=1:runTime
        fprintf('  trial %2d of %2d\n',trial,runTime);

        % Keep A fixed and independently generate xStar, b, and sampling.
        prob=make_suitesparse_problem(matrixFile,delta,...
            caseSeed+trial,true);

        % The access budgets act as safeguards; the displayed epochs are
        % replaced below by the prescribed method-specific definitions.
        cfgVRRK=experiment_config("quick");
        cfgVRRK.sampling=sampling;
        cfgVRRK.maxPasses=maxEpoch;
        cfgVRRK.recordPass=1;
        cfgVRRK.divergenceRSE=1e12;

        cfgAVRRK=cfgVRRK;
        cfgAVRRK.maxPasses=2*maxEpoch;
        cfgAVRRK.recordPass=2;

        cfgREK=cfgVRRK;
        cfgREK.maxPasses=maxEpoch;
        cfgREK.recordPass=1;

        cfgRCD=cfgVRRK;
        cfgRCD.maxPasses=maxEpoch;
        cfgRCD.recordPass=1;

        if stopAtTolerance
            cfgVRRK.tol=tolerance;
            cfgAVRRK.tol=tolerance;
            cfgREK.tol=tolerance;
            cfgRCD.tol=tolerance;
        else
            cfgVRRK.tol=0; %#ok<UNRCH>
            cfgAVRRK.tol=0;
            cfgREK.tol=0;
            cfgRCD.tol=0;
        end

        % Couple the row-index sequence of VRRK-AS and AVRRK-AS.
        rng(caseSeed+100000+trial,'twister');
        hVRRK=run_vrrk_as(prob,parVRRK,cfgVRRK);
        rng(caseSeed+100000+trial,'twister');
        hAVRRK=run_avrrk_as(prob,parAVRRK,cfgAVRRK);

        % Run the established norm-sampled REK and RCD implementations.
        rng(caseSeed+200000+trial,'twister');
        hREK=run_rek(prob,cfgREK);
        rng(caseSeed+300000+trial,'twister');
        hRCD=run_rcd(prob,cfgRCD);

        % Apply the epoch definitions used in the manuscript.
        hVRRK.epoch=hVRRK.iter/prob.m;
        hAVRRK.epoch=hAVRRK.iter/prob.m;
        hREK.epoch=hREK.iter/prob.m+hREK.iter/prob.n;
        hRCD.epoch=hRCD.iter/prob.n;

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
        'H','parVRRK','parAVRRK','tuningLog','parameterTable',...
        'matrixName','delta','maxEpoch','runTime');

    % Add this matrix to the combined terminal summary.
    for methodIndex=1:4
        summaryIndex=summaryIndex+1;
        summaryRows(summaryIndex).matrix=matrixName;
        summaryRows(summaryIndex).m=tuneProb.m;
        summaryRows(summaryIndex).n=tuneProb.n;
        summaryRows(summaryIndex).nnzA=tuneProb.nnz;
        summaryRows(summaryIndex).method=string(methodNames{methodIndex});
        summaryRows(summaryIndex).medianRSE=median(finalRSE(:,methodIndex));
        summaryRows(summaryIndex).q25RSE=...
            sample_percentile(finalRSE(:,methodIndex),25);
        summaryRows(summaryIndex).q75RSE=...
            sample_percentile(finalRSE(:,methodIndex),75);
        summaryRows(summaryIndex).medianEpoch=median(finalEpoch(:,methodIndex));
        summaryRows(summaryIndex).medianTime=median(finalTime(:,methodIndex));
    end

    % Arrange the four trajectory collections in the legend order.
    methodRuns={H.VRRK,H.AVRRK,H.REK,H.RCD};

    % Plot RSE against the method-specific epoch definitions.
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
    title(axEpoch,matrixName,'Interpreter','none');
    legend(axEpoch,epochHandles,methodNames,...
        'Location','best','Interpreter','latex');
    set(axEpoch,'FontSize',fontSize);
    if saveFigures
        exportgraphics(figEpoch,fullfile(caseFolder,...
            matrixName+"_epochs.png"),'Resolution',300);
    end

    % Plot median elapsed time and horizontal time-quantile regions.
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
    title(axTime,matrixName,'Interpreter','none');
    legend(axTime,timeHandles,methodNames,...
        'Location','best','Interpreter','latex');
    set(axTime,'FontSize',fontSize);
    if saveFigures
        exportgraphics(figTime,fullfile(caseFolder,...
            matrixName+"_time.png"),'Resolution',300);
    end
end

%%% save and display the combined summary
summaryTable=struct2table(summaryRows);
writetable(summaryTable,fullfile(outputFolder,...
    'suitesparse_full_rank_summary.csv'));
disp(summaryTable);
fprintf('\nCompleted. Results are stored in %s.\n',outputFolder);

function p = select_suite_row_probability(prob,sampling)
%SELECT_SUITE_ROW_PROBABILITY Return the requested row distribution.
if sampling=="row-norm"
    p=prob.rowP;
else
    p=prob.uniformP;
end
end
