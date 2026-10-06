function h = init_history(maxPasses)
%INIT_HISTORY Preallocate algorithm history.
n=ceil(maxPasses)+2;
h.pass=nan(n,1);
h.epoch=nan(n,1);
h.columnEpoch=nan(n,1);
h.rse=nan(n,1);
h.ne=nan(n,1);
h.iter=nan(n,1);
h.time=nan(n,1);
h.count=0;
h.failed=false;
end
