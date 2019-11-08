function curve2 = ml_closecurve(curve)
%ML_CLOSECURVE Make a curve closed.
%   CURVE2 = ML_CLOSECURVE(CURVE) returns a [closed curve] from the [curve]
%   CURVE.
%   
%   See also

%   20-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

curve2 = curve;

if any(curve(1,:)~=curve(end,:))
    curve2(end+1,:) = curve(1,:);
end