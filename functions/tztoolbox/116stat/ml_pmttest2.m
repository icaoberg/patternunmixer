function [pvalue,ts] = ml_pmttest2(s1,s2)

%ML_PMTTEST2: permutation test
%   PVALUE=ML_PMTTEST2(S1,S2) performs permutation test on two 
%   univariate smaples S1 and S2.
%   PVALUE is the p-value

t1=reshape(s1,1,size(s1,1)*size(s1,2));
t2=reshape(s2,1,size(s2,1)*size(s2,2));

tobs=abs(mean(t1)-mean(t2));
merges=[t1 t2];
L1=length(t1);
L2=length(t2);
Lm=length(merges);

B=50000;
sum=0;

for i = 1:B
    randorder=randperm(Lm);
    %if any(randorder(1:L1)>L1)
    ro1=randorder([1:L1]);
    ro2=randorder([L1+1:Lm]);
    per1=merges(ro1);
    per2=merges(ro2);
    T=abs(mean(per1)-mean(per2));
    sum=sum+(T>=tobs);
    %end
    if  mod(i,10000)==0
        i
    end
end

ts=sum;
pvalue=sum/B;
    
