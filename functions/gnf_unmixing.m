function [alpha,concentration] = ...
    gnf_unmixing(objlearned,data,T1,T2)
% Unmixing to find fraction of all fundamental patterns of a well of images
%
% Created: 12-Apr-2008 T. Peng
% Update: 09-Jun-2010 R.F. Murphy - added object outlier detection
% Update: 11-Jul-2010 R.F. Murphy - added reconstruction error thresholding
%                       and checking for absence of objects in an image
%
% Copyright (C) 2008-2010 Murphy Lab, Carnegie Mellon University
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

%% Initialization

if nargin<4
    T2 = inf;
end
if nargin<3
    T1 = 0.;
end

    zscoredfeat = objlearned.zscoredfeat;
    allfeat = objlearned.allfeat;
    allfluo = objlearned.allfluo;
    imid = objlearned.imid;
    classid = objlearned.classid;
    offset = objlearned.offset;
    scale = objlearned.scale;
    objtype = objlearned.objtype;
    pdf = objlearned.pdf;
    objpost = objlearned.objpost;
    ncluster = objlearned.best_ncluster;

    mixfeat = [];
    mixfluo = [];
    miximid = [];
    h = waitbar(0,'Calculating features for the image(s)...');
    for j = 1:length(data.prot)
        proturl = data.prot{j};
        if ~isempty(data.dna)
            dnaurl = data.dna{j}; 
            [objfeat,objfluo] = gnf_imfeat(proturl,dnaurl,data.protch,data.dnach);
        else
            [objfeat,objfluo] = gnf_imfeatwodna(proturl,data.protch);
        end
        nobjs = length(objfluo);
        mixfeat = [mixfeat;objfeat];
        mixfluo = [mixfluo;objfluo];
        miximid = [miximid;j*ones(nobjs,1)];
        waitbar(j/length(data.prot),h);
    end
    close(h);

    %% Object distribution learning
[objnumdistr,objfluodistr,objnumfrac,objfluofrac] = ...
    gnf_objdistr(objpost,allfluo,imid,classid);
Asum = sum(objfluofrac');

npatterns = size(objnumdistr,1);

% check for no objects found
if size(mixfeat,1)==0
    alpha = cell(3,1);
    for k=1:3
% set all fractions to zero except the last one (the "unknown" pattern)
        alpha{k} = [zeros(npatterns,1); 1]';
    end
    concentration = zeros(npatterns,1);
    return
end

%% object recognition and filtering

zsoremixfeat = ml_zscore(mixfeat,offset,scale);
% Recoginize object types
testobjtype = tz_classify(zsoremixfeat,zscoredfeat,objtype,'dist',1);
testobjdistr = hist(testobjtype,ncluster);
imnum = length(unique(miximid));

% filter objects that are too dissimilar to their type
if T1>0. 
    outlieridx = [];
    for k = 1:ncluster
        clusterfeat = zscoredfeat(objtype==k,:);
        mu = mean(clusterfeat);
        C = cov(clusterfeat);
        isoutlier = tp_isoutlier(zsoremixfeat(testobjtype==k,:),mu,C,T1);
        typekobjidx = find(testobjtype==k);
        outlieridx = [outlieridx;typekobjidx(isoutlier)];
    end
    % Remove outliers
    outlierobjfrac = length(outlieridx)/length(testobjtype);
    outlierfluo = mixfluo(outlieridx);
    outlierfrac = sum(outlierfluo)/sum(mixfluo);
    testobjtype(outlieridx) = [];
    mixfluo(outlieridx) = [];
    imid(outlieridx) = [];
else
    outlierobjfrac = 0.;
    outlierfrac = 0.;
end

%% Unmixing
for k = 1:ncluster
totalfluo(k) = sum(mixfluo(testobjtype==k))/imnum;
end

alpha = cell(3,1);

% Approach 1: Object number frequency linear unmixing
b = testobjdistr/imnum;
A = objnumdistr;
A = gnf_normalize(A,2);
b = gnf_normalize(b,2);
alpha{1} = gnf_lsqlin(A',b')';
rmserr_1 = sqrt(sum((b'-A'*alpha{1}').^2));
% Approach 2: Object number fraction multinomial linear unmixing
alpha{2}(1) = tz_fitmnom2(testobjdistr,pdf(1,:),pdf(2,:));
alpha{2}(2) = 1 - alpha{2}(1);
% Approach 3: Fluorescence fraction concentration estimation
% and unmixing
% normalize before unmixing?
w = gnf_solvemix(objfluofrac',totalfluo');
rmserr_3 = sqrt(sum((totalfluo'-objfluofrac'*w).^2))/sum(totalfluo);
if rmserr_1 > T2
    outlierobjfrac = 1;
end
if rmserr_3 > T2
    outlierfrac = 1;
end
concentration = w;
alpha{3} = gnf_normalize(w.*Asum',1)';

%% renormalize unmixed results to account for outlier fraction
alpha{1} = alpha{1} * (1-outlierobjfrac);
alpha{1}(npatterns+1) = outlierobjfrac;
alpha{2} = alpha{2} * (1-outlierobjfrac);
alpha{2}(npatterns+1) = outlierobjfrac;
alpha{3} = alpha{3} * (1-outlierfrac);
alpha{3}(npatterns+1) = outlierfrac;

