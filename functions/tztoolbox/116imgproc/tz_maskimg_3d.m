function img2=tz_maskimg_3d(img,maskimg)
%TZ_MASKIMG_3D Crop 3D image by an mask.
%   IMG2 = TZ_MASKIMG_3D(IMG,MASKIMG)
%   
%   See also

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

%function img2=tz_maskimg_3d(img,maskimg)
%
%OVERVIEW:
%   

nslice=size(img,3);

for i=1:nslice
    simg=img(:,:,i);
    simg(find(maskimg==0))=0;
    img2(:,:,i)=simg;
end

