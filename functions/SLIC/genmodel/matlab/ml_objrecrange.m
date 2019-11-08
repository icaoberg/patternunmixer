function range=ml_objrecrange(obj)
%ML_OBJRECRANGE Minimal rectangle enclosing an object.
%   RANGE = ML_OBJRECRANGE(OBJ) returns a 1xM vector, of which each column
%   is the expansion of corresponding column of the object OBJ, which is a
%   NxM matrix, where N is the number of pixels and M is the number of
%   dimensions.
%   
%   See also

%   18-Sep-2005 Initial write T. Zhao
%   Copyright (c) Murphy Lab, Carnegie Mellon University

% range(1)=max(obj(:,1))-min(obj(:,1))+1;
% range(2)=max(obj(:,2))-min(obj(:,2))+1;
range=max(obj,[],1)-min(obj,[],1)+1;
