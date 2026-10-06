function h = add_history(h,prob,x,work,iter,elapsed)
%ADD_HISTORY Append accuracy diagnostics.
h.count=h.count+1;
q=h.count;
h.pass(q)=work/prob.nnz;
h.epoch(q)=iter/prob.m;
h.columnEpoch(q)=iter/prob.n;
h.iter(q)=iter;
h.time(q)=elapsed;
h.rse(q)=norm(x-prob.xStar)^2/max(norm(prob.xStar)^2,eps);
h.ne(q)=norm(prob.A'*(prob.A*x-prob.b))/max(norm(prob.A'*prob.b),eps);
end
