function idx = ml_objinimg(obj,imgsize)
%ML_OBJINIMG Test if an object has pixels outside of an image.
%   IDX = ML_OBJINIMG(OBJ,IMGSIZE) returns the indices of points in the
%   [object] or [point array] which are outside of an image with image size
%   IMGSIZE. 
%   
%   See also

%   15-May-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

idx = find(obj(:,1)<1 | obj(:,2)<1 | ...
    obj(:,1)>imgsize(1) | obj(:,2)>imgsize(2));
