function prob = build_problem(A,b,xStar,rStar)
%BUILD_PROBLEM Cache norms and sampling distributions.
A=sparse(A);
[m,n]=size(A);
validateattributes(b,{'numeric'},{'real','finite','column','numel',m});
validateattributes(xStar,{'numeric'},{'real','finite','column','numel',n});
validateattributes(rStar,{'numeric'},{'real','finite','column','numel',m});
assert(isreal(A) && all(isfinite(nonzeros(A))) && nnz(A)>0,...
    'VRRK:InvalidMatrix','A must be real, finite, and nonzero.');
rowNorm2=full(sum(abs(A).^2,2));
colNorm2=full(sum(abs(A).^2,1))';
% Zero rows receive zero probability under row-norm sampling.
if any(rowNorm2==0)
    warning('VRRK:ZeroRows',...
        'The matrix contains zero rows; norm sampling assigns them probability zero.');
end

% Cache all matrix statistics shared by the four solvers.
prob.A=A;
prob.b=b;
prob.xStar=xStar;
prob.rStar=rStar;
prob.m=m;
prob.n=n;
prob.nnz=nnz(A);
prob.rowNorm2=rowNorm2;
prob.colNorm2=colNorm2;
prob.uniformP=ones(m,1)/m;
prob.rowP=rowNorm2/sum(rowNorm2);
prob.colP=colNorm2/sum(colNorm2);
prob.Luniform=max(rowNorm2./prob.uniformP);
end
