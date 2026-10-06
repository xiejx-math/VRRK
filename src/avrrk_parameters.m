function values = avrrk_parameters(prob,par,sampling)
%AVRRK_PARAMETERS Resolve empirical or Theorem 4.1 singleton parameters.
% Theory mode uses a dense SVD and is intended for small validation problems.
if sampling=="row-norm"
    p=prob.rowP;
else
    p=prob.uniformP;
end
assert(all(p>0),'VRRK:ZeroProbability',...
    'Remove zero rows before resolving AVRRK parameters.');
L=max(prob.rowNorm2./p);
mode="empirical";
if isfield(par,'mode')
    mode=string(par.mode);
end
if mode=="theory"
    assert(numel(prob.A)<=5e6,'VRRK:TheorySize',...
        'Theory mode uses a dense SVD; use a small validation matrix.');
    singularValues=svd(full(prob.A),'econ');
    threshold=max(size(prob.A))*eps(max(singularValues));
    positive=singularValues(singularValues>threshold);
    assert(~isempty(positive),'VRRK:ZeroRank','A must have positive numerical rank.');
    mu=min(positive)^2;
    eta=1/(1+min(p));
    % Q=P for singleton sampling, hence lambda_max(Q/P)=1.
    alpha=min((1-eta)^2/(64*L*(2-eta)),1/(4*max(positive)^2));
    theta1=min(0.25,sqrt(alpha*mu/2));
    theta2=0.5;
    gamma=alpha/(2*(theta1+alpha*mu));
    beta=1-gamma*mu;
else
    assert(mode=="empirical",'VRRK:ParameterMode','Unknown AVRRK parameter mode.');
    alpha=par.cAlpha/L;
    theta1=min(0.25,par.cTheta*sqrt(alpha));
    theta2=min(par.theta2,1-theta1);
    gamma=par.cGamma*alpha/2;
    beta=min(1,max(0,1-par.cBeta*gamma));
    eta=par.eta;
    mu=NaN;
end
assert(all(isfinite([alpha theta1 theta2 gamma beta eta])) && ...
    alpha>0 && theta1>0 && theta2>=0 && theta1+theta2<=1 && ...
    gamma>0 && beta>=0 && beta<=1 && eta>0 && eta<=1,...
    'VRRK:InvalidParameters','Invalid realized AVRRK parameters.');
values=struct('mode',mode,'alpha',alpha,'theta1',theta1,...
    'theta2',theta2,'gamma',gamma,'beta',beta,'eta',eta,'mu',mu);
end
