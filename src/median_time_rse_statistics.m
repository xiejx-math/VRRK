function [medianTime,stats,epochGrid] = median_time_rse_statistics(runs)
%MEDIAN_TIME_RSE_STATISTICS Separate time/RSE medians at aligned epochs.
% Reuse the epoch plot's grid and RSE bands. Never zero-pad missing tails.
% Select time and RSE from the same previous recorded checkpoint.
[epochGrid,stats]=curve_statistics(runs,'epoch');
alignedTime=nan(numel(runs),numel(epochGrid));
for trial=1:numel(runs)
    h=runs{trial};
    assert(numel(h.time)==numel(h.epoch) && ...
        all(isfinite(h.time)) && all(diff(h.time)>=0),...
        'VRRK:InvalidHistory','Time must be finite, monotone, and aligned with epochs.');
    [epochs,keep]=unique(h.epoch,'stable');
    times=h.time(keep);
    assert(all(diff(epochs)>0) && epochGrid(1)>=epochs(1) && ...
        epochGrid(end)<=epochs(end),'VRRK:InvalidHistory',...
        'Epoch checkpoints must increase and cover the shared grid.');
    if numel(epochs)==1
        alignedTime(trial,:)=times(1);
    else
        alignedTime(trial,:)=interp1(epochs,times,epochGrid,'previous');
    end
end
medianTime=median(alignedTime,1);
end
