function img = ml_crd2img(pts,param)
%ML_CRD2IMG Convert coordinate shape into an image.
%   IMG = ML_CRD2IMG(PTS) returns an image that contains the [curve] PTS.
%   
%   IMG = ML_CRD2IMG(PTS,PARAM) specifies the parameters of convertion.
%   Currently it contains a field 'tz_obj2img' which has two subfields, 
%   'imgsize' and 'mode'. These two subfileds are parameters for 2nd and
%   3rd arguments for the function ML_OBJ2IMG. 
%   
%   See also

%   31-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('1 or 2 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
end

obj2imgParameters.imgsize = [];
obj2imgParameters.mode = [];
param = ml_initparam(param,struct('tz_obj2img',obj2imgParameters));
pts = ml_showpts_2d(pts,'ln',0);
img = ml_obj2img(pts,param.tz_obj2img.imgsize,param.tz_obj2img.mode);

