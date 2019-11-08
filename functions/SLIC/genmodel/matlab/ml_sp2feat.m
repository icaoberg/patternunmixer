function [knots,coefs] = ml_sp2feat(sp)
%ML_SP2FEAT Convert B-spline structure into features.
%   KNOTS = ML_SP2FEAT(SP) returns the knots of the spline SP.
%   
%   [KNOTS,COEFS] = ML_SP2FEAT(...) also returns coeffiecients. Both KNOTS
%   and COEFS are scalars or row vectors. KNOTS could also be empty.
%   
%   See also

%   28-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

if sp.number-sp.order>0
    knots = sp.knots(sp.order+1:sp.number);
else
    knots = [];
end

coefs = sp.coefs;
