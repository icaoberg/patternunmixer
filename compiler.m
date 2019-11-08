% Compile the ummix GUI

% Author: Tao Peng, Ting Zhao, Ivan E. Cao-Berg and Robert F. Murphy
% Created: April 6, 2009
%
% Copyright (C) 2009 Center for Bioimage Informatics/Murphy Lab
% Carnegie Mellon University
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.
%
% For additional information visit http://murphylab.web.cmu.edu or
% send email to murphy@cmu.edu

addpath( genpath(pwd) );

if ismac
	mkdkir('unmix-macosx-v1.0');
	mcc -m PUnmix -d unmix-macosx-v1.0 -o punmix
	!cp -r data unmix-macosx-v1.0;
	!cp LICENSE unmix-macosx-v1.0;
	!rm punmix-macosx-v1.0/*.c;
	!rm punmix-macosx-v1.0/*.prj;
	!rm punmix-macosx-v1.0/*.log;
	!mv punmix-macosx-v1.0/readme.txt punmix-macosx-v1.0/README
elseif ispc
	mkdir('punmix-win32-v1.0');
    mcc -m PUnmix -d punmix-win32-v1.0 -o punmix
	!copy data punmix-win32-v1.0;
    !copy LICENSE punmix-win32-v1.0;
	!del punmix-win32-v1.0\*.c;
	!del punmix-win32-v1.0\*.prj;
	!del punmix-win32-v1.0\*.log;
	!mv punmix-win32-v1.0\readme.txt punmix-win32-v1.0\README
elseif isunix
	mkdir('punmix-glnx86-v1.0');
    mcc -m PUnmix -d punmix-glnx86-v1.0 -o punmix
	!cp -r data punmix-glnx86-v1.0;
	!cp LICENSE punmix-glnx86-v1.0;
	!rm punmix-glnx86-v1.0/*.c;
	!rm punmix-glnx86-v1.0/*.prj;
	!rm punmix-glnx86-v1.0/*.log;
	!mv punmix-glnx86-v1.0/readme.txt punmix-glnx86-v1.0/README
else
	warning('Unknown or unsupported OS');
end
