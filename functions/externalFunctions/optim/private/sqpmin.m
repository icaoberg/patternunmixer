function [x,val,output,exitflag,lambda] = sqpmin(c,H,xstart,A,b,lb,ub,verb,...
    options,defaultopt,computeLambda,varargin)
%SQPMIN	Solve quadratic problems with box constraints or linear equalities
%
% Locate local solution to
%
%        min { q(x) = .5x'Hx + c'x : l <= x <= u}. 
%
%                            or
%
%       min { q(x) = .5x'Hx + c'x : Ax = b},
%
%
% where H is sparse symmetric matrix. (may be virtual),
%
% x = sqpmin(c,H,xstart,options) return the minimizer of the
% quadratic function q(x) subject to any bounds indicated in
% the named parameter list options. xstart is the starting point.
%
% x = sqpmin(c,H,xstart,options,A,b) solves the linearly
% constrained problem, min { q(x) = .5x'Hx + c'x : Ax = b}.
%
% [x,val] =  sqpmin(c,H,xstart,A,b,ub,lb) returns the value of 
% the quadratic objective function at the solution.
%
% [x,val,gopt] = sqpmin(c,H,xstart,options, ...) returns a measure
% of first-order optimality.
%
% [x,val,gopt,it] = sqpmin(c,H,xstart,options,...) returns
% number of iterations used.
%
% [x,val,gopt,it,npcg] =  sqpmin(c,H,xstart,options,...) returns
% total number of conjugate gradient iterations used.
%
% [x,val,gopt,it,npcg,exitflag] = sqpmin(c,H,xstart,options,...) returns
% termination code.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/10 21:50:52 $

if nargin < 2
  error('optim:sqpmin:NotEnoughInputs','sqpmin requires at least 2 arguments.')
end
if nargin <=2, xstart = []; end
n = length(c); 

if isempty(lb), lb = -inf*ones(n,1); end
if isempty(ub),ub = inf*ones(n,1); end
arg = (ub >= 1e10); arg2 = (lb <= -1e10);
ub(arg) = inf;
lb(arg2) = -inf;
if any(ub == lb) 
   error('optim:sqpmin:EqualLowerUpperBnd', ...
         ['Equal upper and lower bounds not permitted in this large-scale method.\n' ...
          'Use equality constraints and the medium-scale method instead.'])
elseif min(ub-lb) <= 0
   error('optim:sqpmin:InconsistentBnds','Inconsistent bounds.')
end
if isempty(xstart) 
    xstart = startx(ub,lb); 
end
% If components of initial x not within bounds, set those components
% of initial point to a "box-centered" point
xinitOutOfBounds_idx = xstart < lb | xstart > ub;
if any(xinitOutOfBounds_idx)
    xstart = startx(ub,lb,xstart,xinitOutOfBounds_idx);
end

numberOfVariables = n;
if n == 0
  error('optim:sqpmin:InvalidN','n must be positive.')
end

%   INITIALIZATIONS
lambda.lower = [];
lambda.upper = [];
lambda.eqlin = [];  
lambda.ineqlin = [];  % This won't change because no inequalities.
val = []; gopt=[];
output = [];
it = 1; 
if nargin < 5, A = []; end
if isempty(A)
   
   % Box-constrained problem
   
   [x,val,gopt,it,npcg,exitflag,lambda,msg]=sqpbox(c,H,lb,ub,xstart,options,defaultopt,...
       numberOfVariables,verb,computeLambda,varargin{:});
  lambda.ineqlin = []; lambda.eqlin = [];
  output.iterations = it;   
  output.algorithm = 'large-scale: reflective trust-region';
  output.firstorderopt = gopt;
  output.cgiterations = npcg;
  output.message = msg;
else
   if ((max(lb) > -inf) | (min(ub) < inf))
      error('optim:sqpmin:InvalidConstraints', ...
            'sqpmin doesn''t handle both box constraints and Ax = b.');
   else
      % Equality constrained problem
      [mA,nA] = size(A);
      if nargin < 6, b = zeros(mA,1); end
      if isempty(b), b = zeros(mA,1); end
      
      mtxmpy = optimget(options,'HessMult',defaultopt,'fast') ;
      if isempty(mtxmpy)
          mtxmpy = @hmult;
      end
      % Note we pass options in so some values are different for PPCGR than for SQPBOX
      [x,po,npcg,pgnrm,exitflag,lambda,msg]=ppcgr(c,H,mtxmpy,A,b,options,defaultopt,verb,computeLambda,...,
          [],[],[],[],[],[],varargin{:});
      if exitflag == -10
         % ppcgr aborted
         return
      end
      
      it = 1; % number of iterations reported in output.cgiterations
      % Calculate final objective value 
      w = feval(mtxmpy,H,x,varargin{:}); 
      gopt = pgnrm;
      val = x'*(c + .5*w);
      
      output.iterations = it; 
      output.algorithm = 'large-scale: projective preconditioned conjugate gradients';
      output.firstorderopt = gopt;
      output.cgiterations = npcg;
      output.message = msg;
   end
end





