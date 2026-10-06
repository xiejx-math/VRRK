function value = sample_percentile(data,percentage,dimension)
%SAMPLE_PERCENTILE Midpoint-rank percentile (MATLAB prctile convention).
% Order statistics have probabilities (i-0.5)/N, with clipped endpoints.
% NaNs are rejected; Inf is allowed for unattained time-to-accuracy targets.
if nargin<3
    dimension=find(size(data)~=1,1);
    if isempty(dimension)
        dimension=1;
    end
end
validateattributes(data,{'numeric'},{'real','nonempty'});
assert(~any(isnan(data(:))),'VRRK:NaNPercentile','NaN samples are not supported.');
validateattributes(percentage,{'numeric'},{'scalar','finite','>=',0,'<=',100});
validateattributes(dimension,{'numeric'},{'scalar','integer','positive','<=',ndims(data)});
order=1:ndims(data);
order([1 dimension])=order([dimension 1]);
sorted=sort(permute(data,order),1);
count=size(sorted,1);
position=min(count,max(1,count*percentage/100+0.5));
lower=floor(position);
upper=ceil(position);
weight=position-lower;
subs=repmat({':'},1,ndims(sorted));
subs{1}=lower;
a=sorted(subs{:});
subs{1}=upper;
b=sorted(subs{:});
if lower==upper
    value=a;
else
    value=(1-weight)*a+weight*b;
end
value=ipermute(value,order);
end
