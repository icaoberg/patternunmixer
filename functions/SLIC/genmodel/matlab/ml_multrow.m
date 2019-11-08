function B = ml_multrow(A,x)
%ML_MULTROW Multiply the a row vector to each row of a matrix.
%   ML_MULTROW(A,X) returns a matrix in which each row is the array product
%   between the corresponding row in matrix A and the row vector X. A and X
%   must have the same number of columns if X is not a scalar.
%   
%   See also ML_ADDROW

%   13-Mar-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('Exactly 2 arguments are required')
end

if length(x(:))==1
    B = A*x;
else
    if(size(A,2)~=size(x,2))
        error('The matrix and the vector must have the same number of columns');
    end
    B=A.*(ones(size(A,1),1)*x);
end