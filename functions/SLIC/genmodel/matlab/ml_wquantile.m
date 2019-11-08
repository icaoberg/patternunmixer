function y = ml_wquantile(x,p)
%ML_WQUANTILE Weighted quantile of a sample.
%   Y = ML_WQUANTILE(X,P) returns the P quantile of the weights in X, which 
%   is a row vector or column vector.
%   
%   See also

%   26-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

if size(x,1)>1
    x = x';
end

[s,idx] = sort(x,2,'descend');

cs = cumsum(s);
cs = cs/cs(end)>p;
[tmp,sidx] = max(cs);

y = s(sidx);






