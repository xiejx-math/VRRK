function [commonGrid,stats] = curve_statistics(runs,xField)
%CURVE_STATISTICS Compute pointwise statistics after trajectory alignment.
assert(iscell(runs) && ~isempty(runs),'VRRK:EmptyRuns','At least one trial is required.');
assert(ismember(string(xField),["epoch","time"]),'VRRK:InvalidAxis',...
    'The horizontal coordinate must be epoch or time.');
if string(xField)=="epoch"
    getX=@(h) h.epoch;
else
    getX=@(h) h.time;
end

% Use only the horizontal interval covered by every trial.
xEnd=inf;
xStart=-inf;
sameGrid=true;
firstGrid=getX(runs{1});
for ii=1:numel(runs)
    x=getX(runs{ii});
    y=runs{ii}.rse;
    assert(isvector(x) && ~isempty(x) && numel(x)==numel(y) && ...
        all(isfinite(x)) && all(diff(x)>=0),'VRRK:InvalidHistory',...
        'Coordinates must be finite, nondecreasing, and aligned with RSE.');
    xStart=max(xStart,x(1));
    xEnd=min(xEnd,x(end));
    sameGrid=sameGrid && isequal(x(:),firstGrid(:));
end
assert(xStart<=xEnd,'VRRK:NoCommonRange','Trials have no shared horizontal interval.');
if sameGrid
    % Preserve every original checkpoint when trials already align.
    commonGrid=unique(firstGrid(:)','stable');
else
    commonGrid=linspace(xStart,xEnd,501);
end

% Align recorded RSE with a previous-value step function.
alignedRSE=nan(numel(runs),numel(commonGrid));
for ii=1:numel(runs)
    x=getX(runs{ii});
    y=runs{ii}.rse;
    assert(all(isfinite(y)) && all(y>=0),'VRRK:InvalidHistory',...
        'Nonfinite or negative RSE must be investigated before plotting.');
    [x,keep]=unique(x,'stable');
    y=y(keep);
    if numel(x)==1
        alignedRSE(ii,:)=y(1);
    else
        alignedRSE(ii,:)=interp1(x,y,commonGrid,'previous');
    end
end

% Compute statistics in the original RSE scale, then plot on a log axis.
stats.minimum=min(alignedRSE,[],1);
stats.q25=sample_percentile(alignedRSE,25,1);
stats.median=median(alignedRSE,1);
stats.q75=sample_percentile(alignedRSE,75,1);
stats.maximum=max(alignedRSE,[],1);
end
