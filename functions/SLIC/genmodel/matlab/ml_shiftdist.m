function dists2 = ml_shiftdist(dists,angle)
%ML_SHIFTDIST Shift the distance of hit points to a centain angle.
%   DISTS2 = ML_SHIFTDIST(DISTS,ANGLE) returns a vector of hit point
%   distances that are shifted from DISTS by ANGLE. This function assumes
%   that DISTS is started from angle 0 and the step is one degree.
%   
%   See also

%   10-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

dists2 = circshift(dists,[0 -round(angle)]);