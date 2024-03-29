function xIndices = classifyBoundsOnVars(lb,ub,nVar)
%classifyBoundsOnVars Helper function that identifies variables
% that are fixed, and that have finite bounds. 
 
%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2007/11/09 20:56:42 $

% Set empty vector bounds to vectors of +Inf or -Inf
if isempty(lb)
    lb = -Inf(nVar,1);
end
if isempty(ub)
    ub = Inf(nVar,1);
end
% Check for NaN
if any(isnan(lb)) || any(isnan(ub))
    error('optim:classifyBoundsOnVars:NaNBounds', ...
        'Bounds lb and ub must not be NaN.')
end
% Check for +Inf lb and -Inf ub
if any(lb == Inf)
    error('optim:classifyBoundsOnVars:PlusInfLb', ...
        'Lower bounds lb cannot be +Inf.')
end
if any(ub == -Inf)
    error('optim:classifyBoundsOnVars:MinusInfUb', ...
        'Upper bounds ub cannot be -Inf.')
end
% Check for fixed variables equal to Inf or -Inf
if any(lb(:) == ub(:) & isinf(lb(:)))
    error('optim:classifyBoundsOnVars:FixedVarsAtInf', ...
        'Bounds lb and ub cannot be equal and infinite.')
end

% Fixed variables
xIndices.fixed = equalFloat(lb,ub,eps);

% Finite lower and upper bounds; exclude fixed variables
xIndices.finiteLb = ~xIndices.fixed & isfinite(lb(:));
xIndices.finiteUb = ~xIndices.fixed & isfinite(ub(:));

%-------------------------------------------------------------------------
function isEqual_idx = equalFloat(v1,v2,tolerance)
%equalFloat Helper function that compares two vectors
% using a relative difference and returns a boolean
% vector.

% Indices for which both v1 and v2 are finite
finiteRange_idx = isfinite(v1(:)) & isfinite(v2(:));

% Indices at which v1 and v2 are (i) finite and (ii) equal in a 
% floating point sense
isEqual_idx = finiteRange_idx & ...
    abs(v1(:)-v2(:)) <= tolerance*max( 1,max(abs(v1(:)),abs(v2(:))) );
