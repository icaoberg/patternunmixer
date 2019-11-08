function [objnumdistr,objfluodistr,fluo_concen_k,avgfluofrac]...
    = gnf_objdistr(objpost,fluo,imid,classid)
%GNF_OBJDISTR calculates distribution of objects and object fluorescence
%fraction of fundamental patterns
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

for class = 1:max(classid)
        imnum = length(unique(imid(classid==class)));
        objnum(class,:) = sum(objpost(classid==class,:));
        objnumdistr(class,:) = objnum(class,:)/imnum;
        fluopost = objpost.*repmat(fluo,[1,size(objpost,2)]);
        objfluo(class,:) = sum(fluopost(classid==class,:));
        objfluodistr(class,:) = objfluo(class,:)/imnum;
end

objfluodistr(:,1) = [];

mitoconcen = [412 242 142 84 49 29]';
lysoconcen = [214 153 109 78 56 40]';

mitok = sum(mitoconcen.*objfluodistr(1,:))./sum(mitoconcen.^2);
lysok = sum(lysoconcen.*objfluodistr(2,:))./sum(lysoconcen.^2);

%mitoarray = repmat(mitoconcen,[1,size(objfluodistr,2)]);
%lysoarray = repmat(lysoconcen,[1,size(objfluodistr,2)]);

%mitok = sum(mitoarray.*objfluodistr(1,:))./sum(mitoarray.^2);
%lysok = sum(lysoarray.*objfluodistr(2,:))./sum(lysoarray.^2);
fluo_concen_k = [mitok' lysok'];

%% regression with average fluorescence fraction and a single coefficient
objfluofrac = gnf_normalize(objfluodistr,2);
avgfluofrac = squeeze(mean(objfluofrac));
fluosum = squeeze(sum(sum(objfluodistr)));
coefsum(1) = sum(sum(mitoconcen*avgfluofrac(:,1)'));
coefsum(2) = sum(sum(lysoconcen*avgfluofrac(:,2)'));
fluo_concen_k = fluosum./coefsum';
