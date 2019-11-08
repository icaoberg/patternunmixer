function gnf_plotobjfractions(objlearned)
% Plot fractions of object types in each fundamental pattern of a learned
% model
%
% Created: 10-Jul-2010 R.F. Murphy
%
% Copyright (C) 2010 Murphy Lab, Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty ofpatt
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

%% Initialization

    zscoredfeat = objlearned.zscoredfeat;
    allfeat = objlearned.allfeat;
    allfluo = objlearned.allfluo;
    imid = objlearned.imid;
    classid = objlearned.classid;
    offset = objlearned.offset;
    scale = objlearned.scale;
    objtype = objlearned.objtype;
    pdf = objlearned.pdf;
    objpost = objlearned.objpost;
    ncluster = objlearned.best_ncluster;

%% Object distribution learning
[objnumdistr,objfluodistr,objnumfrac,objfluofrac] = ...
    gnf_objdistr(objpost,allfluo,imid,classid);
Asum = sum(objfluofrac');
figure(101)
bar(objnumfrac')
xlabel('Object type')
ylabel('Fraction of objects')
legend(objlearned.patternlist)
figure(102)
bar(objfluofrac')
xlabel('Object type')
ylabel('Fraction of fluorescence')
legend(objlearned.patternlist)

