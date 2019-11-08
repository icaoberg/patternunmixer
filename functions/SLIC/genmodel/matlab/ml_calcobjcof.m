function cof=ml_calcobjcof(obj)
%ML_CALCOBJCOF Calculate COF of an object.
%   COF = ML_CALCOBJCOF(OBJ) calculate COF of an object. The object is 
%   a matrix with 2 columns or 3 columns. The first two columns are X and
%   Y coordinates. The third column is the vector of gray levels, which
%   are all ones if the third column does not exist.

%   27-JUN-2004 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

if nargin < 1
    error('Exactly 1 argument is required')
end

if size(obj,2)==2
    cof(1)=mean(obj(:,1));
    cof(2)=mean(obj(:,2));   
else
    cof(1)=sum(obj(:,1).*obj(:,3))/sum(obj(:,3));
    cof(2)=sum(obj(:,2).*obj(:,3))/sum(obj(:,3));
end
