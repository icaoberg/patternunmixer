function alpha=gnf_solvemix(A,b,norm)
%TZ_SOLVEMIX Solve linear equation by pseudo-inverse.
%   ALPHA = TZ_SOLVEMIX(A,b)
%   
%   See also

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

%function alpha=tz_solvemix(y,p)
%
%OVERVIEW:
%   solve linear equations approximately
%PARAMETERS:
%   y - mixture pattern
%   p - fundamental patterns
%RETURN:
%   alpha - coeffiecients
%DESCRIPTION:
%   equation: Ax=y
%   solution: x=inv(T(A)A)T(A)y
%HISTORY:
%   07-DEC-2004 Initial write TINGZ
%

if ~exist('norm','var')
    norm = 'nonormalize';
end

alpha = zeros(size(A,2),1);
tmpalpha = inv(A'*A)*A'*b;
allindex = 1:length(alpha);

while any(tmpalpha<0)
    [minalpha,minind] = min(tmpalpha);
    A(:,minind)=[];
    allindex(minind)=[];
    if isempty(allindex)
        warning('solution failed');
        return;
    end
    tmpalpha = inv(A'*A)*A'*b;
end

alpha(allindex)=tmpalpha;
if strcmp(norm,'normalize')
    alpha=alpha/sum(alpha);
end
