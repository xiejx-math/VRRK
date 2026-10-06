function [rseGrid,stats] = time_curve_statistics(runs,rseGrid)
%TIME_CURVE_STATISTICS First RECORDED attainment times, without interpolation.
% Unreached user-specified targets are represented by Inf, never extrapolated.
% By default, use only targets reached by every trial of this method.
assert(~isempty(runs),'VRRK:EmptyRuns','At least one trial is required.');
if nargin<2
    initial=cellfun(@(h) h.rse(1),runs);
    terminal=cellfun(@(h) min(h.rse),runs);
    assert(all(isfinite([initial(:);terminal(:)])),...
        'VRRK:InvalidHistory','Nonfinite RSE must be investigated before plotting.');
    high=min(initial);
    low=max(terminal);
    assert(high>=low && high>0,'VRRK:NoCommonRange',...
        'Trials do not share a positive attained RSE interval.');
    low=max(low,realmin);
    rseGrid=logspace(log10(high),log10(low),501);
    rseGrid(1)=high;
    rseGrid(end)=low;
end
validateattributes(rseGrid,{'numeric'},{'vector','positive','finite'});
rseGrid=rseGrid(:)';
times=inf(numel(runs),numel(rseGrid));
for trial=1:numel(runs)
    h=runs{trial};
    assert(numel(h.time)==numel(h.rse) && all(diff(h.time)>=0),...
        'VRRK:InvalidHistory','Time and RSE histories must align.');
    for target=1:numel(rseGrid)
        first=find(h.rse<=rseGrid(target),1,'first');
        if ~isempty(first)
            times(trial,target)=h.time(first);
        end
    end
end
stats.minimum=min(times,[],1);
stats.q25=sample_percentile(times,25,1);
stats.median=median(times,1);
stats.q75=sample_percentile(times,75,1);
stats.maximum=max(times,[],1);
stats.reachedCount=sum(isfinite(times),1);
stats.timeAtRSE=times;
end
