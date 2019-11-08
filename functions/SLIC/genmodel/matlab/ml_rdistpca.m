function [f,avgratio] = ml_rdistpca(combcellcodes,param)
%ML_RDISTPCA Train the the PCA model for distance ratios.
%   F = ML_RDISTPCA(COMBCELLCODES) returns the trained PCA model for the
%   input cell array of cell codes COMBCELLCODES. F is a structure with the
%   following fields:
%       'startangle' - see below.
%       'anglestep' - see below.
%       'minratio' - minimal ratios. It will be ignored if it is empty. If it
%          is 0, the minimal ratio will be found automatically.
%       'stat' - a pdf.
%   
%   F = ML_RDISTPCA(COMBCELLCODES,PARAM) specifies parameters for training.
%   PARAM is a structure that has the following fields:
%       'startangle' - a string that determines how to align the ratio
%           vector. 'cell' means the major angle of the cell and 'nuc'
%           means the major angle of the nucleus.
%       'anglestep' - step of the angles. It must be an integer.
%       'minratio' - if it less than 1, the ratio will be obtained from the
%           data.
%       'ml_estpdf' - parameters for the function ML_ESTPDF
%
%   [F,AVGRATIO] = ML_RDISTPCA(...) also returns the average shape ratio.
%
%   See also

%   10-Jan-2006 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 1
    error('Exactly 1 argument is required')
end

if ~exist('param','var')
    param = struct([]);
end

if ~isfield(param,'anglestep')
    if length(combcellcodes)>360
        param.anglestep=1;
    else
        param.anglestep = round(359/length(combcellcodes)+0.5);
    end
end

param = ml_initparam(param, ...
    struct('startangle','nuc','minratio',[],'ml_estpdf',struct([])));

for i=1:length(combcellcodes)
    cellcode = combcellcodes{i};
    rdist = cellcode.nucelldist./cellcode.nucdist;

    switch param.startangle
        case 'cell'
            startAngle = cellcode.cellmangle;
        case 'nuc'
            startAngle = cellcode.nucmangle;
    end
    rdist = rdist((0:param.anglestep:359)+1);
    normrdist(i,:) = ml_shiftdist(rdist,startAngle);
end

f.startangle = param.startangle;
f.anglestep = param.anglestep;
if ~isempty(param.minratio)
    if param.minratio<1
        f.minratio = min(normrdist(:));
    else
        f.minratio = param.minratio;
    end
else
    model.minratio = [];
end

avgratio = mean(normrdist,1);


%normrdist = ml_addrow(normrdist,-avgratio);
%f.avgratio = avgratio;

% if isfield(param.ml_estpdf.transform.param,'ncomp')
%     param.ml_estpdf.mu = zeros(1,param.ml_estpdf.transform.param.ncomp);
% else
%     param.ml_estpdf.mu = zeros(1,size(normrdist,2));
% end

f.stat = ml_estpdf(normrdist,param.ml_estpdf);

