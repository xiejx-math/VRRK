function h = run_avrrk_as(prob,par,cfg,indices)
%RUN_AVRRK_AS Run accelerated VRRK-AS with practical parameter coupling.

% Unpack the problem and select the row-sampling distribution.
A=prob.A;
b=prob.b;
n=prob.n;
if cfg.sampling=="row-norm"
    p=prob.rowP;
else
    p=prob.uniformP;
end

% The default is empirical tuning, not the theoretical parameter schedule.
values=avrrk_parameters(prob,par,cfg.sampling);
alpha=values.alpha;
theta1=values.theta1;
theta2=values.theta2;
gamma=values.gamma;
beta=values.beta;
eta=values.eta;

% Cache the sampling distribution before timing starts.
cdf=cumsum(p);
cdf(end)=1;
if nargin>=4
    validateattributes(indices,{'numeric'},{'vector','integer','positive'});
    assert(all(indices<=prob.m) && all(p(indices)>0),...
        'VRRK:InvalidIndices','Supplied indices must have positive probability.');
end

% Initialize primal, momentum, and residual-tracking variables.
x=zeros(n,1);
y=x;
yprev=x;
z=x;
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

% Perform accelerated randomized row updates.
while work<cfg.maxPasses*prob.nnz && k<maxIterations
    k=k+1;
    if nargin>=4 && k<=numel(indices)
        i=indices(k);
    else
        i=draw_cdf(cdf);
    end

    % Mix momentum states and form the variance-reduced direction.
    ai=A(i,:)';
    xmix=theta1*z+theta2*yprev+(1-theta1-theta2)*y;
    dg=full(ai'*xmix)-b(i)-rbar(i);
    g=e+(dg/p(i))*ai;

    % Update the main and momentum sequences.
    ynew=xmix-alpha*g;
    znew=beta*z+(1-beta)*xmix+(gamma/alpha)*(ynew-xmix);

    % Refresh the residual tracker using the second row inner product.
    dr=full(ai'*yprev)-b(i)-rbar(i);
    e=e+eta*dr*ai;
    rbar(i)=rbar(i)+eta*dr;
    yprev=y;
    y=ynew;
    z=znew;
    x=y;
    work=work+2*nnz(ai);
% Two row inner products are evaluated.

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
h.reportedSequence="y";
h.tolerance=cfg.tol;
h.samplingIncluded=true;
if nargin>=4 && ~isempty(indices)
    h.samplingIncluded=false;
end


% Store realized parameters and trim unused history allocation.
h.parameterMode=values.mode;
h.alpha=alpha;
h.theta1=theta1;
h.theta2=theta2;
h.gamma=gamma;
h.beta=beta;
h.eta=eta;
h=finalize_history(h);
end
