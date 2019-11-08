function r = cubici1(f2,f1,c2,c1,dx)
%CUBICI1 Cubicly interpolates 2 points and gradients to estimate minimum.
%
%   This function uses cubic interpolation and the values of two 
%   points and their gradients in order to estimate the minimum of a 
%   a function along a line.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/06/14 05:18:57 $

if isinf(f2), f2 = 1/eps; end
z = 3*(f1-f2)/dx+c1+c2;
w = real(sqrt(z*z-c1*c2));
r = dx*((z+w-c1)/(c2-c1+2*w));
