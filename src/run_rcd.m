function h = run_rcd(prob,cfg)
%RUN_RCD Randomized coordinate descent; use only for full column rank.
A=prob.A;
b=prob.b;
x=zeros(prob.n,1);
residual=-b;
% Cache distributions; random draws themselves are included in elapsed time.
colCDF=cumsum(prob.colP);
colCDF(end)=1;
work=0;
k=0;
next=cfg.recordPass;
h=init_history(cfg.maxPasses);
% Initial diagnostics and initialization are outside the timed iteration.
h=add_history(h,prob,x,work,k,0);
h.stopReason="budget";
maxIterations=inf;
if isfield(cfg,'maxIterations')
    maxIterations=cfg.maxIterations;
end
validateattributes(cfg.recordPass,{'numeric'},{'scalar','positive','finite'});
validateattributes(cfg.maxPasses,{'numeric'},{'scalar','positive','finite'});
validateattributes(maxIterations,{'numeric'},{'scalar','positive'});
excludedTime=0;
timer=tic;
while work<cfg.maxPasses*prob.nnz && k<maxIterations
    k=k+1;
    j=draw_cdf(colCDF);
    aj=A(:,j);
    step=full(aj'*residual)/prob.colNorm2(j);
    x(j)=x(j)-step;
    residual=residual-step*aj;
    work=work+nnz(aj);
    if work/prob.nnz>=next
        [h,excludedTime,stop]=record_solver_state(...
            h,prob,x,work,k,timer,excludedTime,cfg);
        next=work/prob.nnz+cfg.recordPass;
        if stop
            break;
        end
    end
end

% Always record the actual final iterate, even between diagnostic checkpoints.
if h.iter(h.count)~=k
    h=record_solver_state(h,prob,x,work,k,timer,excludedTime,cfg);
end
h.xFinal=x;
h.tolerance=cfg.tol;
h.samplingIncluded=true;

h=finalize_history(h);
h.epoch=h.iter/prob.n;
end
