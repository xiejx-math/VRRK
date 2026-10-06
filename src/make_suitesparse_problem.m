function prob = make_suitesparse_problem(file,delta,seed,fullColumnRank)
%MAKE_SUITESPARSE_PROBLEM Read Matrix Market and create inconsistent b.
A=mmread_simple(file);
% Structurally zero rows contain no information about x and cannot be used
% by norm-proportional row sampling. Removing them leaves the least-squares
% solution unchanged and avoids trivial, permanently inconsistent equations.
nonzeroRows=full(sum(spones(A),2))>0;
removedZeroRows=sum(~nonzeroRows);
A=A(nonzeroRows,:);
[m,n]=size(A);
rng(seed,'twister');
if fullColumnRank
    xStar=randn(n,1);
else
    xStar=A'*randn(m,1);
end
xStar=xStar/max(norm(xStar),eps);
z=randn(m,1);
y=lsqminnorm(A,z);
be=z-A*y;
be=delta*max(norm(A*xStar),eps)*be/max(norm(be),eps);
b=A*xStar+be;
prob=build_problem(A,b,xStar,-be);
prob.removedZeroRows=removedZeroRows;
end
