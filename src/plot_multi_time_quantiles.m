function handles = plot_multi_time_quantiles(ax,runs,lineColors,...
    lightColors,darkColors,lineWidth)
%PLOT_MULTI_TIME_QUANTILES Median time versus median RSE at aligned epochs.
% Bands describe RSE variation, not time variation. Both plots use identical
% RSE statistics; only their horizontal coordinates differ.
hold(ax,'on');
set(ax,'YScale','log');

methodNumber=numel(runs);
timeGrids=cell(methodNumber,1);
statistics=cell(methodNumber,1);

% Align checkpoints by epoch and take time and RSE medians separately.
for ii=1:methodNumber
    [timeGrids{ii},statistics{ii}]=median_time_rse_statistics(runs{ii});
end

% Draw every minimum-to-maximum RSE region first.
for ii=1:methodNumber
    t=timeGrids{ii};
    stats=statistics{ii};
    fill(ax,[t fliplr(t)],...
        [stats.minimum fliplr(stats.maximum)],...
        lightColors(ii,:),...
        'EdgeColor','none',...
        'FaceAlpha',0.35,...
        'HandleVisibility','off');
end

% Draw every 25th-to-75th percentile RSE region next.
for ii=1:methodNumber
    t=timeGrids{ii};
    stats=statistics{ii};
    fill(ax,[t fliplr(t)],...
        [stats.q25 fliplr(stats.q75)],...
        darkColors(ii,:),...
        'EdgeColor','none',...
        'FaceAlpha',0.55,...
        'HandleVisibility','off');
end

% Draw all median RSE curves above every shaded region.
handles=gobjects(methodNumber,1);
for ii=1:methodNumber
    handles(ii)=plot(ax,timeGrids{ii},statistics{ii}.median,...
        'Color',lineColors(ii,:),...
        'LineWidth',lineWidth);
end
end
