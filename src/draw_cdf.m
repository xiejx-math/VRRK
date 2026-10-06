function index = draw_cdf(cdf)
%DRAW_CDF Draw one sample using a cached CDF and binary search.
% Strict comparison prevents selection of zero-probability entries.
u=rand;
left=1;
right=numel(cdf);
while left<right
    middle=floor((left+right)/2);
    if u<cdf(middle)
        right=middle;
    else
        left=middle+1;
    end
end
index=left;
end
