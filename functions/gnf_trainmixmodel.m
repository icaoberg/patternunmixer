function [modelFile, retrain] = gnf_trainmixmodel(data,name,modelpath,maxObjects)
%Train mixture model on collections of images of fundamental patterns

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

if nargin < 2 || isempty(name)
    name = 'New_Model';
end

% Saving directories
datapath = './data/train';
if ~exist(modelpath,'dir')
    mkdir(modelpath)
end

% Calculate features
% featfile = strcat(modelpath,'/features.mat');
% if exist(featfile,'file')
%     load(featfile)
% else
    npatterns = length(data);
    allfeat = [];
    allfluo = [];
    patternid = [];
    imgid = [];
    for i = 1:npatterns
        h = waitbar(0,strcat('Calculating features for ''',data(i).name,''' image(s)'));
        for j = 1:length(data(i).prot)
            proturl = data(i).prot{j};
            if ~isempty(data(i).dna)
                dnaurl = data(i).dna{j};
                [objfeat,objfluo] = gnf_imfeat(proturl,dnaurl,data(i).protch,data(i).dnach);
            else
                [objfeat,objfluo] = gnf_imfeatwodna(proturl,data(i).protch);
            end
            nobjs = length(objfluo);
            allfeat = [allfeat;objfeat];
            allfluo = [allfluo;objfluo];
            imgid = [imgid;j*ones(nobjs,1)];
            patternid = [patternid;i*ones(nobjs,1)];
            waitbar(j/length(data(i).prot),h);
        end
        close(h);
    end
%     save(featfile,'allfeat','allfluo','imgid','patternid');
% end

% Learn object types
[model, retrain] = gnf_objectlearn(allfeat,allfluo,imgid,patternid,maxObjects);
model.name = name;
model.patternlist = cell(length(data),1);
for i = 1:length(data)
    model.patternlist{i} = data(i).name;
end
modelFile = strcat(modelpath,filesep,'model_',regexprep(name,'\s*|:|\\|/|\.|\''|"','_'),'.mat');
save(modelFile,'model');
