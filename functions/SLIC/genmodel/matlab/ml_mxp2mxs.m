function shape2 = ml_mxp2mxs(shape)
%ML_MXP2MXS Conver medial axis spline into medial axis.
%   SHAPE2 = ML_MXP2MXS(SHAPE) returns the medial axis representation of
%   the shape which is originally represented by spline shape SHAPE.
%   
%   See also

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

x=(0:shape.length-1)/(shape.length-1);

shape2.format = 'mxs';
shape2.medaxis = [(1:shape.length)',spval(shape.spmedaxis,x)'];
shape2.width = spval(shape.spwidth,x);

