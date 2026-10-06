function prob = make_synthetic_problem(m,n,r,kappa,delta,seed,solutionType)
%MAKE_SYNTHETIC_PROBLEM Controlled-rank inconsistent least-squares system.
if nargin<7
    solutionType="random";
end
validateattributes(m,{'numeric'},{'scalar','integer','positive'});
validateattributes(n,{'numeric'},{'scalar','integer','positive'});
validateattributes(r,{'numeric'},{'scalar','integer','positive','<=',min(m,n)});
validateattributes(kappa,{'numeric'},{'scalar','finite','>=',1});
validateattributes(delta,{'numeric'},{'scalar','finite','positive'});
assert(r<m,'VRRK:NoResidualSpace',...
    'An inconsistent system requires rank(A)<m.');
assert(ismember(string(solutionType),["random","tail"]),...
    'VRRK:SolutionType','solutionType must be random or tail.');
rng(seed,'twister');
[U,~]=qr(randn(m,r),0);
[V,~]=qr(randn(n,r),0);
if r==1
    s=1;
else
    s=logspace(0,-log10(kappa),r)';
end
A=U*(s.*V');
if solutionType=="tail"
    q=min(5,r);
    xStar=V(:,r-q+1:r)*(ones(q,1)/sqrt(q));
else
    xStar=V*randn(r,1);
    xStar=xStar/norm(xStar);
end
z=randn(m,1);
be=z-U*(U'*z);
be=delta*max(norm(A*xStar),eps)*be/max(norm(be),eps);
b=A*xStar+be;
prob=build_problem(sparse(A),b,xStar,-be);
prob.rank=r;
prob.kappa=kappa;
prob.delta=delta;
prob.solutionType=string(solutionType);
end

