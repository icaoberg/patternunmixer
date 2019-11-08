function [objfeat,objfluo] = gnf_imfeatwodna(proturl,protch)
%GNF_IMFEAT calculate object level features within an image without DNA
%related features

% Copyright (C) 2009 Center for Bioimage Informatics/Murphy Lab
% Carnegie Mellon University
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

% Check if feature already exist
featfile = strcat(proturl,'.mat');
if exist(featfile,'file')
    load(strcat(proturl,'.mat'));
    return
end

% Load images and preprocess
protimg = imread(proturl);

s = size(protimg);

if(length(s) > 2 && s(3) >= protch)
    protimg = protimg(:,:,protch);
end
procimg = ml_preprocess(double(protimg),[],'ml','yesbgsub');
if isempty(procimg) 
    objfeat = [];
    objfluo = [];
    return
end

% Calculate protein object morphological features
protobjs = ml_findobjs(procimg);

for k = 1:length(protobjs)
    objfeat(k,:) = tz_objfeat...
        (protobjs{k},[],struct('featset',{{'mor','skl'}}));
    objfluo(k) = sum(protobjs{k}(:,end));
end
objfluo = objfluo';
objfeat(:,2:3) = [];
