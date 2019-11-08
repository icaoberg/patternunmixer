function img = ml_mxs2img(medaxis,width)
%ML_MXS2IMG Convert medial axis representation into an image.
%   IMG = ML_MXS2IMG(MEDAXIS,WIDTH) returns an image that contains the
%   shape represented by medial axis MEDAXIS and width WIDTH.
%   
%   See also

%   01-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

img = ml_crd2img(ml_mxs2crd(medaxis,width));