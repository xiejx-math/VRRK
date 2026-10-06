function score = history_score(h)
%HISTORY_SCORE Area under log-error curve; invalid candidates score Inf.
if h.failed || any(~isfinite(h.rse)) || isempty(h.rse)
    score=inf;
    return;
end
x=h.pass/max(h.pass(end),eps);
y=min(log10(max(h.rse,realmin)),0);
score=trapz(x,y)+0.25*log10(max(max(h.rse),1));
end

