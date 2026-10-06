function outputFolders = Replot_saved_results(resultFolder,visible)
%REPLOT_SAVED_RESULTS Redraw saved trials without rerunning any algorithm.
% Replot_saved_results opens a folder picker. Select one experiment folder
% (containing kappa/rank/matrix subfolders), or a single case folder.
% Replot_saved_results(resultFolder,false) exports without showing figures.
% Time curves: median(time) versus median(RSE) at shared epoch checkpoints.
% Both plots use identical RSE bands. Missing tails are never zero-padded.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'));
if nargin<1 || strlength(string(resultFolder))==0
    resultFolder=uigetdir(fullfile(root,'results'),'Select saved experiment folder');
    if isequal(resultFolder,0)
        outputFolders={};
        return;
    end
end
if nargin<2
    visible=true;
end
assert(isfolder(resultFolder),'The selected result folder does not exist.');

% Collect saved histories, including cases below the selected folder.
files=dir(fullfile(resultFolder,'all_trial_histories.mat'));
nested=dir(fullfile(resultFolder,'**','all_trial_histories.mat'));
files=[files; nested];
paths=arrayfun(@(f) fullfile(f.folder,f.name),files,'UniformOutput',false);
[~,keep]=unique(paths,'stable');
files=files(keep);
assert(~isempty(files),'No all_trial_histories.mat files were found.');
outputFolders=cell(numel(files),1);

% Edit these display settings to customize the figures.
lineWidth=2;
fontSize=14;
colors=[0 .447 .741; .85 .325 .098; .466 .674 .188; .494 .184 .556];
keys={'VRRK','AVRRK','REK','RCD'};
labels={'VRRK-AS','AVRRK-AS','REK','RCD'};
visibility='off';
if visible
    visibility='on';
end
for caseIndex=1:numel(files)
    sourceFile=fullfile(files(caseIndex).folder,files(caseIndex).name);
    data=load(sourceFile,'H');
    assert(isfield(data,'H'),'Saved file does not contain H: %s',sourceFile);
    present=cellfun(@(key) isfield(data.H,key),keys);
    caseKeys=keys(present);
    methodNames=labels(present);
    assert(~isempty(caseKeys),'No supported method histories in %s',sourceFile);
    methodRuns=cellfun(@(key) data.H.(key),caseKeys,'UniformOutput',false);
    lineColors=colors(present,:);
    lightColors=.8+.2*lineColors;
    darkColors=.4+.6*lineColors;

    % Verify that changing the horizontal coordinate preserves all RSE data.
    plotStatistics=struct();
    for methodIndex=1:numel(caseKeys)
        runs=methodRuns{methodIndex};
        [medianTime,stats,epoch]=median_time_rse_statistics(runs);
        [epochCheck,statsCheck]=curve_statistics(runs,'epoch');
        assert(isequal(epoch,epochCheck) && isequal(stats,statsCheck));
        plotStatistics.(caseKeys{methodIndex})=struct('epoch',epoch,...
            'medianTime',medianTime,'rseStatistics',stats,'numTrials',numel(runs));
    end

    % New folders preserve the original figures and raw numerical results.
    out=new_result_folder(files(caseIndex).folder,'replot_double_median');
    outputFolders{caseIndex}=out;
    save(fullfile(out,'plot_statistics.mat'),'plotStatistics','sourceFile');
    [~,caseTitle]=fileparts(files(caseIndex).folder);
    for plotIndex=1:2
        fig=figure('Color','w','Visible',visibility);
        ax=axes(fig);
        hold(ax,'on');
        set(ax,'YScale','log');
        if plotIndex==1
            handles=plot_multi_epoch_quantiles(ax,methodRuns,...
                lineColors,lightColors,darkColors,lineWidth);
            xLabel='Epochs';
            stem='RSE_epochs';
        else
            handles=plot_multi_time_quantiles(ax,methodRuns,...
                lineColors,lightColors,darkColors,lineWidth);
            xLabel='Elapsed time (seconds)';
            stem='RSE_elapsed_time';
        end
        grid(ax,'on');
        box(ax,'on');
        xlabel(ax,xLabel);
        ylabel(ax,'RSE');
        title(ax,caseTitle,'Interpreter','none');
        legend(ax,handles,methodNames,'Location','best');
        set(ax,'FontSize',fontSize);
        ax.Toolbar.Visible='off';
        drawnow;
        exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',300);
        savefig(fig,fullfile(out,[stem '.fig']));
        if ~visible
            close(fig);
        end
    end
    fprintf('Replotted %s\nSaved to %s\n',sourceFile,out);
end
end

