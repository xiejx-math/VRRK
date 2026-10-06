function handles = plot_multi_epoch_quantiles(ax,runs,lineColors,...
    lightColors,darkColors,lineWidth)
%PLOT_MULTI_EPOCH_QUANTILES Plot uncertainty bands for several methods.

methodNumber=numel(runs);
grids=cell(methodNumber,1);
statistics=cell(methodNumber,1);

% Align the independent runs of each method on a common epoch grid.
for ii=1:methodNumber
    [grids{ii},statistics{ii}]=curve_statistics(runs{ii},'epoch');
end

% Draw every minimum-to-maximum region before drawing any median line.
for ii=1:methodNumber
    x=grids{ii};
    stats=statistics{ii};
    fill(ax,[x fliplr(x)],...
        [stats.minimum fliplr(stats.maximum)],...
        lightColors(ii,:),...
        'EdgeColor','none',...
        'FaceAlpha',0.35,...
        'HandleVisibility','off');
end

% Draw every interquartile region over the corresponding outer region.
for ii=1:methodNumber
    x=grids{ii};
    stats=statistics{ii};
    fill(ax,[x fliplr(x)],...
        [stats.q25 fliplr(stats.q75)],...
        darkColors(ii,:),...
        'EdgeColor','none',...
        'FaceAlpha',0.55,...
        'HandleVisibility','off');
end

% Draw all median curves last so that no shaded region hides a curve.
handles=gobjects(methodNumber,1);
for ii=1:methodNumber
    handles(ii)=plot(ax,grids{ii},statistics{ii}.median,...
        'Color',lineColors(ii,:),...
        'LineWidth',lineWidth);
end
end
