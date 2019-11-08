function [objfeat,objfluo,cellnum] = gnf_imfeat(prot,dna,protch,dnach,imsize)
% Calculate features of objects in single protein image
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

[dnalabeled,num_dna] = bwlabel(dna>0);
for k = 1:num_dna
    if bwarea(dnalabeled==k) < 400
        dna(dnalabeled==k) = 0;
    end
end
[dnalabeled,num_dna] = bwlabel(dna>0);
cellnum = num_dna;

if cellnum <= 5
    objfeat = [];
    objfluo = [];
    return
end

for k = 1:num_dna
    dnaimg = dna;
    dnaimg(dnalabeled~=k) = 0;
    dna_m00 = ml_imgmoments(dnaimg,0,0) ;
    dna_m01 = ml_imgmoments(dnaimg,0,1) ;
    dna_m10 = ml_imgmoments(dnaimg,1,0) ;
    dnacof(k,:) = [dna_m01/dna_m00,dna_m10/dna_m00] ;
end

% Calculate protein object morphological features
protobjs = ml_findobjs(prot);

for k = 1:length(protobjs)
    objfeat(k,:) = tz_objfeat...
        (protobjs{k},[],struct('featset',{{'mor','skl'}}));
    objfluo(k) = sum(protobjs{k}(:,end));
    protcof(k,:) = ml_calcobjcof(protobjs{k});
end
objfluo = objfluo';

% Calculate two DNA related features
[PROTX,DNAX] = meshgrid(protcof(:,1),dnacof(:,1));
[PROTY,DNAY] = meshgrid(protcof(:,2),dnacof(:,2));
distmap = sqrt((PROTX-DNAX).^2+(PROTY-DNAY).^2);
[distfeat,protdnaid] = min(distmap);

for k = 1:length(protobjs)
    protimg = ml_obj2img(protobjs{k},imsize);
%     objsize = bwarea(protimg>0);
    objsize = objfeat(k,1);
    dnamsk = dnalabeled==protdnaid(k);
    overlapfeat(k) = bwarea((protimg&dnamsk)>0)/objsize;
end

dnafeat = [distfeat',overlapfeat'];
objfeat(:,2:3) = dnafeat;