function obj = ml_gaussobj(sigma)
%ML_GAUSSOBJ An object from Gaussian distribution.
%   OBJ = ML_GAUSSOBJ(SIGMA) returns an object that is extracted from a
%   2D Gaussian distribution which has covariance matrix SIGMA. The object
%   contains no less than 95% energy of the Gaussian distribution.
%   
%   See also

%   26-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

img = ml_gaussimg(sigma);

y = ml_wquantile(img(:),0.95);
img(img<y) = 0;

imageSize = size(img);

[r,c]=find(img>0);
obj = [r,c,img(sub2ind(imageSize,r,c))];