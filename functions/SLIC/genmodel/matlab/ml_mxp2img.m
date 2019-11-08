function img = ml_mxp2img(shape)
%ML_MXP2IMG Convert medial axis spline shape into an image.
%   IMG = ML_MXP2IMG(SHAPE) returns an image that contains the medial axis
%   spline shape SHAPE.
%   
%   See also TZ_IMG2MXP

%   01-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

shape2 = ml_mxp2mxs(shape);
img = ml_mxs2img(shape2.medaxis,shape2.width);