function pts2 = ml_deformpts_2d(pts,A)
%ML_DEFORMPTS_2D Transform 2D points.
%   PTS2 = ML_DEFORMPTS_2D(PTS,A) returns an array of 2d points which is
%   the tranformation of PTS according to transformation matrix A.
%   
%   See also

%   07-Apr-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

npt=size(pts,1);

order=-(size(A,1)-1)/2;

qpts=ml_expandpts(pts,order);

pts2=qpts*A;
pts2=pts2(:,1:2);
