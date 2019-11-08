function img2=ml_mainobjcontour(img)
%ML_MAINOBJCONTOUR Find the boundary of the biggest object in an image.
%   IMG2 = ML_MAINOBJCONTOUR(IMG) returns an image that contains the
%   boundary of the biggest object in the image IMG.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if isempty(img)
    img2=[];
    return;
end

img=double(img>0);

img2=imfill(img,'holes');
imgedge=bwperim(img2,4);

img2=double(img2)-double(imgedge);
obj=ml_findmainobj_bw(img2,4);

img2=ml_obj2img(obj,size(img));
img2=img2==0;
img2=[ones(1,size(img2,2));img2;ones(1,size(img2,2))];
img2=[ones(size(img2,1),1),img2,ones(size(img2,1),1)];
img2=bwperim(img2,4);
img2(:,1)=[];
img2(:,end)=[];
img2(1,:)=[];
img2(end,:)=[];
