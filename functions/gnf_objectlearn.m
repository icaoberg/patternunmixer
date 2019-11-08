function [objectlearned, fewMaxClusters]  = gnf_objectlearn(allfeat,allfluo,imid,classid,MaxCluster)
% GNF_OBJECTLEARN

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

%% Clustering
[zscoredfeat,offset,scale] = ml_zscore(allfeat);
if(isempty(MaxCluster))
    MaxCluster = 35;
end
objtypeall = cell(MaxCluster,1);
pdfall = cell(MaxCluster,1);
objpostall = cell(MaxCluster,1);
h = waitbar(0,'Learning mixture model...');
for ncluster = 2:MaxCluster
    [objpostall{ncluster},pdfall{ncluster}] = gnf_objcluster...
        (zscoredfeat,imid,classid,'kmeans',...
        ncluster,'nosave','./data/results/train');
    [objtypeall{ncluster},non] = find(objpostall{ncluster}');
    [aic(ncluster),bic(ncluster)] = tz_aicbic(zscoredfeat,objtypeall{ncluster},'euc');
    waitbar(ncluster/MaxCluster,h);
end
close(h)
[aic_min,best_ncluster] = min(aic(2:end));
best_ncluster = best_ncluster + 1;
fewMaxClusters = 0;
if(MaxCluster - best_ncluster <= 4)
    fewMaxClusters = 1;
end
objtype = objtypeall{best_ncluster};
objpost = objpostall{best_ncluster};
pdf = pdfall{best_ncluster};

objectlearned = struct(...
    'allfeat',allfeat,...
    'zscoredfeat',zscoredfeat,...
    'allfluo',allfluo,...
    'imid',imid,...
    'classid',classid,...
    'aic',aic,...
    'bic',bic,...
    'offset',offset,...
    'scale',scale,...
    'objtype',objtype,...
    'pdf',pdf,...
    'objpost',objpost,...
    'best_ncluster',best_ncluster);
try
    rmdir('./data/results','s')
catch exception
    disp('Matlab Error in rmdir: Could not delete the temp directory:');
    disp('./data/results');
end
