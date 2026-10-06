function h = run_vrrk_as(prob,par,cfg,indices)
%RUN_VRRK_AS Run VRRK-AS with singleton row sampling.

% Unpack the problem and select the row-sampling distribution.
A=prob.A;
b=prob.b;
n=prob.n;
if cfg.sampling=="row-norm"
    p=prob.rowP;
else
    p=prob.uniformP;
end

% Scale the dimensionless tuning parameter by the sampling smoothness.
L=max(prob.rowNorm2./max(p,eps));
alpha=par.cAlpha/L;

% Cache the sampling distribution before timing starts.
cdf=cumsum(p);
cdf(end)=1;
if nargin>=4
    validateattributes(indices,{'numeric'},{'vector','integer','positive'});
    assert(all(indices<=prob.m) && all(p(indices)>0),...
        'VRRK:InvalidIndices','Supplied indices must have positive probability.');
end

% Initialize the iterate, residual tracker, and aggregated direction.
x=zeros(n,1);
rbar=-b;
e=A'*rbar;
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

% Perform randomized row updates until convergence or budget exhaustion.
while work<cfg.maxPasses*prob.nnz && k<maxIterations
    k=k+1;
    if nargin>=4 && k<=numel(indices)
        i=indices(k);
    else
        i=draw_cdf(cdf);
    end

    % Form the variance-reduced direction and update its control variate.
    ai=A(i,:)';
    d=full(ai'*x)-b(i)-rbar(i);
    g=e+(d/p(i))*ai;
    x=x-alpha*g;
    e=e+d*ai;
    rbar(i)=rbar(i)+d;
    work=work+nnz(ai);

    % Record diagnostics approximately once per matrix-access pass.
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
if nargin>=4 && ~isempty(indices)
    h.samplingIncluded=false;
end


% Store the realized stepsize and remove unused history allocation.
h.alpha=alpha;
h=finalize_history(h);
end
