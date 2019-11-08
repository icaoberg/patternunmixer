function [aic,bic] = tz_aicbic(x,label,distfun)
%TZ_AICBIC Calculate aic and bic values of clusters.
%   AIC = TZ_AICBIC(X,LABEL,DISTFUN) returns the AIC score of the clusters
%   of X. LABEL is a vector of the cluster labels. DISTFUN is the distance
%   function: 'euc' for Euclidean distance and 'mah' for Mahalonobis
%   distance.
%   
%   [AIC,BIC] = TZ_AICBIC(...) also returns BIC values.
%   
%   See also

%   24-Dec-2005 Initial write T. Zhao
%   Copyright (c) Center for Bioimage Informatics, CMU

if nargin < 3
    error('Exactly 3 arguments are required')
end

clusternum = ml_combfeats2mcf([1:length(label)]',label);

switch distfun
    case 'euc'
        [aic,bic] = ml_aicbic_euc(x,clusternum);
    case 'mah'
        [aic,bic] = ml_aicbic_mah(x,clusternum);
    case 'euc2'
        ncluster = length(clusternum);
        feat = x;
        removeidx = [];
        for i=1:ncluster
            idx = clusternum{i};
            mu = mean(feat(idx,:),1);
            feat(idx,:) = ml_addrow(feat(idx,:),-mu);
            if length(idx)==1
                removeidx = [removeidx;idx];
            end
        end
        feat(removeidx,:) = [];
        poolvar = var(feat(:),1);
        x = x/sqrt(poolvar);
        [aic,bic] = ml_aicbic_euc(x,clusternum);
        addlk = length(x(:))*log(poolvar)/2;
        aic = aic+addlk;
        bic = bic+addlk;       
    otherwise
        error('Unrecognized distance function');
end
