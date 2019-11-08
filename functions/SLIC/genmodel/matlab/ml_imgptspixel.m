function ps=ml_imgptspixel(img,pts)
%ML_IMGPTSPIXEL Get pixel values at specified points.
%   PS = ML_IMGPTSPIXEL(IMG,PTS) returns a vector of pixel values from
%   the image IMG. PS(I) is the gray level of the pixel at position
%   [PTS(I,1),PTS(I,2)] in IMG.

%   ??-???-???? Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 2
    error('Exactly 2 arguments are required')
end

imgsize=size(img);
pts(pts(:,1)<=0 | pts(:,1)>imgsize(1),:)=[];
pts(pts(:,2)<=0 | pts(:,2)>imgsize(2),:)=[];

if isempty(pts)
    error('empty points');
end

imgsize=size(img);

idx=sub2ind(imgsize,pts(:,1),pts(:,2));

ps=img(idx);
