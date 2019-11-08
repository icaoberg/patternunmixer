function pts2 = ml_moldshape(pts,dists)
%ML_MOLDSHAPE Relocate points according to the distance difference.
%   PTS2 = ML_MOLDSHAPE(PTS,DISTS) returns an array of points that are form
%   the points PTS, which are relocated according to the distance
%   differences DISTS. PTS and DISTS must have the same number of rows.
%   
%   See also

%   11-Apr-2005 Initial write  T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end 

[tmp,s]=min(abs(pts(pts(:,1)>0,2)));

if s>1
    tpts=[pts(s:end,:);pts(1:s-1,:)];
else
    tpts=pts;
end
orgdists=sqrt(sum(tpts.^2,2));
rdist = (orgdists+dists)./orgdists;
pts2(:,1)=tpts(:,1).*rdist;
pts2(:,2)=tpts(:,2).*rdist;

% [th,r]=cart2pol(pts(:,1),pts(:,2));
% r = r+dists;
% [pts2(:,1),pts2(:,2)] = pol2cart(th,r);