function idx = sample_discrete(p,count)
%SAMPLE_DISCRETE Toolbox-free sampling with one cached CDF per call.
validateattributes(p,{'numeric'},{'vector','real','finite','nonnegative','nonempty'});
validateattributes(count,{'numeric'},{'scalar','integer','nonnegative','finite'});
assert(sum(p)>0,'VRRK:InvalidProbability','Probabilities must have positive sum.');
cdf=cumsum(p(:)/sum(p));
cdf(end)=1;
idx=zeros(count,1);
for sampleIndex=1:count
    idx(sampleIndex)=draw_cdf(cdf);
end
end

