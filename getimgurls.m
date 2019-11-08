function [filelist,pathname] = getimgurls(dlgtitle,startpath)
% GETIMGURLS Get full path of images selected by user

% Author: Tao Peng (tpeng@cmu.edu)
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

if nargin < 2 || isempty(startpath)
    startpath = pwd;
end

if startpath ~= filesep
    startpath(end+1) = filesep;
end

[filename,pathname] = uigetfile(...
{   '*.bmp','Bitmap Files (*.BMP)';...
    '*.tif','TIFF (*.TIF;*.TIFF)'; ...
    '*.jpg;*.jpeg;*.jpe;*.jfif','JPEG (*.JPG;*.JPEG;*.JPE;*.JFIF)';...
    '*.gif','GIF (*.GIF)';...
    '*.png','PNG (*.PNG)';...
    '*.*','All Files (*.*)'},...
    dlgtitle, startpath, 'MultiSelect','on');

if ~iscell(filename)
    if ~ischar(filename)
        filelist = [];
        return
    end
    filelist{1} = [pathname filename];
    return
else
    filelist = cell(length(filename),1);
    for i = 1:length(filename)
        filelist{i} = [pathname filename{i}];
    end
end
