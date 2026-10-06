function h = finalize_history(h)
%FINALIZE_HISTORY Remove unused preallocation.
f={'pass','epoch','columnEpoch','rse','ne','iter','time'};
for fieldIndex=1:numel(f)
    fieldName=f{fieldIndex};
    h.(fieldName)=h.(fieldName)(1:h.count);
end
end
