function sp = ml_feat2sp(knots,coefs)
%ML_FEAT2SP Convert knots and coeffents into B-spline
%   SP = ML_FEAT2SP(KNOTS,COEFS) returns a B-spline structure with internal
%   nodes KNOTS and coefficients COEFS.
%   
%   See also

%   28-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

sp.form = 'B-';
sp.coefs = coefs;
sp.number = length(coefs);
sp.order = length(coefs)-length(knots);
sp.dim = 1;
sp.knots = [zeros(1,sp.order) knots ones(1,sp.order)];
