function isoutlier = tp_isoutlier(X,mu,C,alpha)
%TP_ISOUTLIER determines outlier of data
%
% Author: T. Peng
% Created: 15-Apr-2009
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

[n,p] = size(X);
Cinv = inv(C);

% Calculate Mahalanobis distances
md2 = zeros(n,1);
for i = 1:n
    xi = X(i,:) - mu;
    md2(i) = xi * Cinv * xi';
end

% Perform Chi-square test
thr = chi2inv(1-alpha,p);
K = length(thr);
isoutlier = false(n,K);
for k = 1:K
    isoutlier(:,k) = md2 > thr(k);
end