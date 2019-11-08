function W = hmult(Hinfo,Y,varargin)
%HMULT	Hessian-matrix product
%
% W = HMULT(Hinfo,Y) An example of a Hessian-matrix product function
% file, e.g. Hinfo is the actual Hessian and so W = Hinfo*Y.
%
% Note: varargin is not used but must be provided in case 
% the objective function has additional problem dependent
% parameters (which will be passed to this routine as well).

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/08/20 16:37:36 $

W = Hinfo*Y;





