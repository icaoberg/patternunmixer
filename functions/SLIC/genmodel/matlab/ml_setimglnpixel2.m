function img2=ml_setimglnpixel2(img,s,a,len)
%ML_SETIMGLNPIXEL2 Obsolete. See ML_SETIMGLNPIXEL2.
%   IMG2 = ML_SETIMGLNPIXEL2(IMG,S,A,LEN) draw a line from the starting
%   [point] S with an angle A and length LEN.
%   
%   See also ML_SETIMGLNPIXEL

%   17-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University


pts=ml_getlinept2(s,a,len);
img2=ml_setimgptspixel(img,pts);
