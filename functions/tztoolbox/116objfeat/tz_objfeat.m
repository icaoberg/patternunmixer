function [feats,names]=tz_objfeat(obj,dnaproc,param)
%TZ_OBJFEAT Calculate 2D object features.
%   FEATS = TZ_OBJFEAT(OBJ,DNAPROC) returns a feature vector of the object
%   OBJ, which is a 3-column matrix, in which the 1st and 2nd column are x 
%   and y coordinates, and the 3rd column contains the gray levels. DNAPROC
%   is an image of DNA channel. In this image, all pixels with intensities
%   above 0 will be considered as DNA and other pixles are background. If 
%   DNAPROC is empty, all DNA-related features will be NaN. 
%   Here is the list of the features
%   (Reference: Object Type Recognition for Automated Analysis of Protein 
%   Subcellular Location. T. Zhao et al. IEEE Trans. Image Proc. 
%   14:1351-1359):
%   1. Number of pixels in object
%   2. Distance between object center of fluorescence(COF) and DNA COF
%   3. Fraction of object pixels overlapping with DNA
%   4. A measure of eccentricity of the object.
%   5. Euler number of the object
%   6. A measure of roundness of the object.
%   7. The length of the object's skeleton.
%   8. The ratio of skeleton length to the area of the convex hull of the
%       skeleton.
%   9. The fraction of object pixels contained within the skeleton.
%   10. The fraction of object fluorescence contained within the skeleton.
%   11. The ratio of the number of branch points in skeleton to length of
%       skeleton.
%   12-24. 13 Haralick texture features (see ML_TEXTURE).
%   25-29. 5 edge features (see ML_IMGEDGEFEATURES).
%   
%   FEATS = TZ_OBJFEAT(OBJ,DNAPROC,PARAM) allows customizing parameters for
%   feature calculation. PARAM is a structure and it has the following
%   fields:
%       'featset' - a cell array that contains the names of feature sets 
%           to calculate. 'mor' corresponds to feature 1 to 6, 'skl' 
%           corresponds to feaure 7-11, 'har' corresponds to feature 12-24 
%           and 'edg' corresponds to feature 25-29.    
%       * 'org_pixsize' - original pixel size (default 0.23um)
%       * 'har_pixsize' - pixel size for haralick texture features 
%            (default 1.15um)
%       * 'har_intbins' - number of gray level (default 256);
%
%       * for calculating haralick texture features only.
%
%   The order of the features will be changed according to PARAM.featset.
%   One way to make sure that you are using the right features is to return
%   another vavariale NAMES (see below).
%
%   [FEATS,NAMES] = TZ_OBJFEAT(...) returns the names of the features in
%   the cell array NAMES: 1.'obj_size', 2.'obj_dna_dist', 
%   3.'obj_dna_overlap', 4.'obj_eccentricity', 5.'obj_euler', 
%   6.'obj_roundness', 7.'obj_skel_len', 8.'obj_skel_hull_area_ratio',
%   9.'obj_skel_obj_area_ratio', 10.'obj_skel_obj_fluor_ratio',
%   11.'obj_skel_branch_per_len', 12.'angular_second_moment', 
%   13.'contrast', 14.'correlation', 15. 'sum_of_squares', 
%   16.'inverse_diff_moment', 17.'sum_avg', 18.'sum_var', 19.'sum_entropy',
%   20.'entropy', 21.'diff_var', 22.'diff_entropy', 
%   23.'info_measure_corr_1', 24.'info_measure_corr_2', 
%   25.'edges:area_fraction', 26.'edges:homogeneity', 
%   27.'edges:direction_maxmin_ratio', 
%   28.'edges:direction_maxnextmax_ratio', 29.'edges:direction_difference'.          
%
%   See also

%   11-NOV-2004 Initial write  T. Zhao
%   30-JUL-2007 Modify T. Peng
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 2
    error('2 or 3 arguments are required')
end

if ~exist('param','var')
    param = struct([]);
    allFeatset = {'mor','skl','har','edg'};
    param = ml_initparam(param,struct('featset',{allFeatset}));
end

objimg = ml_obj2img(obj,[]);
feats = [];
names = {};

for i=1:length(param.featset)
    switch param.featset{i}
        case 'mor'
            objmask = im2bw(objimg) ;
            dnamask   = im2bw(dnaproc) ;
            objsize = bwarea(objmask) ;
            objcof=ml_calcobjcof(obj);
            
            if ~isempty(dnaproc)
                %Calculate DNA COF
                dnaproc_m00 = ml_imgmoments(dnaproc,0,0) ;
                dnaproc_m01 = ml_imgmoments(dnaproc,0,1) ;
                dnaproc_m10 = ml_imgmoments(dnaproc,1,0) ;
                dnacof = [dnaproc_m01/dnaproc_m00 ...
                    dnaproc_m10/dnaproc_m00] ;

                %distance between objcof and dnacof
                obj_dna_dist=sqrt(sum((objcof-dnacof).^2));

                %Fractional overlap with the DNA image
                obj_dna_overlap = bwarea(im2bw( ...
                    ml_obj2img(obj,size(dnaproc))) & dnaproc)/objsize;
            else
                obj_dna_dist = NaN;
                obj_dna_overlap = NaN;
            end
            
            obj_eccent = NaN;
            obj_euler = NaN;
            
            obj_features = ...
                regionprops(uint8(objmask),'Eccentricity','EulerNumber');
            if(~isempty(obj_features))
                obj_eccent = obj_features.Eccentricity;
                obj_euler = obj_features.EulerNumber;
            end
            % Roundness of the object
            objimageperim = bwarea(bwperim(objmask)) ;
            obj_round = (objimageperim^2)/(4*pi*objsize) ;
            feats = [feats objsize obj_dna_dist obj_dna_overlap ...
                obj_eccent obj_euler obj_round] ;
            names = {names{:} 'obj_size' 'obj_dna_dist' ...
                'obj_dna_overlap' 'obj_eccentricity' 'obj_euler' ...
                'obj_roundness'};
        case 'skl'
            if(~isempty(find(objimg)))
                [subfeats, subnames] = ml_objskelfeats( objimg);
            else
                subnames = {'obj_skel_len' ...
                'obj_skel_hull_area_ratio' ...
                'obj_skel_obj_area_ratio' ...
                'obj_skel_obj_fluor_ratio' ...
                'obj_skel_branch_per_len'};
                subfeats = nan(1,length(subnames));
            end
            feats = [feats subfeats];
            names = {names{:} subnames{:}};
        case 'har'
            param = ml_initparam(param,...
                struct('org_pixsize',0.23, ...
                'har_pixsize',1.15,'har_intbins',256));
            
            %tz- 13-Oct-2006
%             [subnames, subfeats] = ...
%                 ml_features(objimg, [], ones(size(objimg)), {'har'});
            %tz--
            
            %tz- 07-Apr-2007
            %tz+ 13-Oct-2006
%             [subnames, subfeats] = ...
%                 ml_features(objimg, [], ones(size(objimg)), ...
%                 {'har'}, 0.23,0,[],256,0.23);
            %tz++
            %tz--
            [subnames, subfeats] = ...
                ml_features(objimg, [], ones(size(objimg)), ...
                {'har'}, param.org_pixsize,[],[],param.har_pixsize, ...
                param.har_intbins);
            %tz+
            
            %tz++
            feats = [feats subfeats];
            names = {names{:} subnames{:}};
        case 'edg'
            [subnames,subfeats] = ml_imgedgefeatures(objimg);
            feats = [feats subfeats];
            names = {names{:} subnames{:}};    
    end
end

% [names, feats, feat_slf] = ...
%     ml_features(objimg, [], ones(size(objimg)), {'har'});


% feats=ml_texture(uint8(objimg/max(objimg(:))*255));
% feats=double(feats);
% feats(isnan(feats))=0;

% [names, feats, slfnames] = ml_imgfeatures(objimg, dnaproc);

%[names, values] = mv_imgobjfeatures(objimg, dnaproc);

% names = {} ;
% values = [] ;
% 
% objmask = im2bw(objimg) ;
% dnamask   = im2bw(dnaproc) ;
% 
% objsize = bwarea(objmask) ;
% objcof=tz_calcobjcof(obj);
% % 
% % Calculate DNA COF
% %
% dnaproc_m00 = mb_imgmoments(dnaproc,0,0) ;
% dnaproc_m01 = mb_imgmoments(dnaproc,0,1) ;
% dnaproc_m10 = mb_imgmoments(dnaproc,1,0) ;
% dnacof = [dnaproc_m10/dnaproc_m00 ...
%                    dnaproc_m01/dnaproc_m00] ;
% 
% %distance between objcof and dnacof
% obj_dna_dist=sqrt(sum((objcof-dnacof).^2));
% 
% 
% %
% % Fractional overlap with the DNA image
% %
% obj_dna_overlap = bwarea(objmask & dnaproc)/objsize;
% 
% %
% % From imfeature()
% %
% obj_features = imfeature(objmask,'Eccentricity','EulerNumber');
% obj_eccent = obj_features.Eccentricity;
% obj_euler = obj_features.EulerNumber;
% 
% %
% % Roundness of the object
% %
% objimageperim = bwarea(bwperim(objmask)) ;
% obj_round = (objimageperim^2)/(4*pi*objsize) ;
% 	
% [edgenames, edgevalues, slfnames] = ml_imgedgefeatures(objimg);
% 
% objlabel=bwlabel(objmask);
% obj_sepnum=max(objlabel(:));
% 
% names = [cellstr('obj_size') ...
%          cellstr('obj_dna_dist') ...
%          cellstr('obj_dna_overlap') ...
%          cellstr('obj_eccentricity') ...
%          cellstr('obj_euler') ...
%          cellstr('obj_roundness') ...
%          cellstr('obj_sepnum') ...
%          edgenames] ;
% feats = [objsize obj_dna_dist obj_dna_overlap obj_eccent obj_euler ...
%   obj_round obj_sepnum edgevalues] ;
