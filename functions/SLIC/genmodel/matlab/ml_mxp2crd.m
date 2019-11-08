function pts = ml_mxp2crd(shape)
%ML_MXP2CRD Converts medial axis spline into coordinates.
%   PTS = ML_MXP2CRD(SHAPE) returns an array of points from the medial axis
%   spline shape SHAPE.
%   
%   See also

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

shape2 = ml_mxp2mxs(shape);
pts = ml_mxs2crd(shape2.medaxis,shape2.width);
