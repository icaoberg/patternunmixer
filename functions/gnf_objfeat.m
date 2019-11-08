function [objfeat,objfluo,cellnum] = gnf_objfeat(path,well,idx)
% Calculate object features of a single protein image
%
% Author: T. Peng
% Created: 12-Apr-2008
% Last Update: 23-Dec-2009
%
% Copyright (C) 2008-2010 Murphy Lab, Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

[image,factor] = gnf_getimage(path,well,idx);
image = double(image)*factor;
[imagedna,factordna] = gnf_getimage(path,well,idx+1);
imagedna = double(imagedna)*factordna;
imagesize = size(image);

TotalFluoThresh = 1.5e7;
if sum(imagedna(:)) < TotalFluoThresh
    objfeat = [];
    objfluo = [];
    cellnum = 0;
    return
end

[procimg,mask,res] = ml_preprocess(image,[],'ml','yesbgsub');
[procdna,mask,res] = ml_preprocess(imagedna,[],'ml','yesbgsub');

[dnas,dnanum] = bwlabel(mask);
if dnanum <= 4
    impurity = imdilate(dnas>0,strel('disk',15));
    imagedna(impurity) = 0;
    [procdna,mask,res] = ml_preprocess(imagedna,[],'ml','yesbgsub');
end

if ~isempty(procimg) && ~isempty(procdna)
    [objfeat,objfluo,cellnum] = gnf_imfeat(procimg,procdna,imagesize);
else
    objfeat = [];
    objfluo = [];
    cellnum = 0;
end