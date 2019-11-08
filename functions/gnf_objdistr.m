function [objnumdistr,objfluodistr,objnumfrac,objfluofrac]...
    = gnf_objdistr(objpost,fluo,cellid,classid)
%GNF_OBJDISTR calculates distribution of objects and object fluorescence
%fraction of fundamental patterns

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

%% Calculate object frequency and fluorescence distribution
nclusters = size(objpost,2);
nclasses = max(classid);
for c = 1:nclasses
    ncells = length(unique(cellid(classid==c)));
    objnum(c,:) = sum(objpost(classid==c,:));
    objnumdistr(c,:) = objnum(c,:)/ncells;
    fluopost = objpost.*repmat(fluo,[1,nclusters]);
    objfluo(c,:) = sum(fluopost(classid==c,:));
    objfluodistr(c,:) = objfluo(c,:)/ncells;
end

%% Calculate object frequency and fluorescence fraction distribution
objnumfrac = gnf_normalize(objnumdistr,2);
objfluofrac = gnf_normalize(objfluodistr,2);
