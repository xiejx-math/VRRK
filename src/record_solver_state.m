function [h,excludedTime,stop] = record_solver_state(h,prob,x,work,k,timer,excludedTime,cfg)
%RECORD_SOLVER_STATE Record diagnostics without charging their cost to solvers.
elapsed=toc(timer)-excludedTime;
diagnosticTimer=tic;
h=add_history(h,prob,x,work,k,elapsed);
value=h.rse(h.count);
stop=false;
if ~isfinite(value) || value>cfg.divergenceRSE
    h.failed=true;
    h.stopReason="diverged";
    stop=true;
elseif cfg.tol>0 && value<=cfg.tol
    h.stopReason="tolerance";
    stop=true;
end
excludedTime=excludedTime+toc(diagnosticTimer);
end
