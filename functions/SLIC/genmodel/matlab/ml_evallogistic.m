function p=ml_evallogistic(x,para)
%ML_EVALLOGISTIC Probability of a logistic model.
%   P = TZ_EVALLOGISTIC(X,PARA) returns the probability of a logitic model with
%   the parameters PARA, in which each column is a set of parameters. X must 
%   be a [feature matrix] and the number of columns must be the same as the 
%   number of rows of PARA.
%   
%   See also

%   ??-OCT-2004 Initial write TINGZ
%   31-OCT-2004 Modified TINGZ
%       - add comments
%   Copyright (c) Murphy Lab, Carnegie Mellon University


y=x*para;

p=1./(1+exp(-y));
