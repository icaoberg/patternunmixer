function [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqnonlin(FUN,x,LB,UB,options,varargin)
%LSQNONLIN solves non-linear least squares problems.
%   LSQNONLIN attempts to solve problems of the form:
%   min  sum {FUN(X).^2}    where X and the values returned by FUN can be   
%    X                      vectors or matrices.
%
%   X = LSQNONLIN(FUN,X0) starts at the matrix X0 and finds a minimum X to 
%   the sum of squares of the functions in FUN. FUN accepts input X 
%   and returns a vector (or matrix) of function values F evaluated
%   at X. NOTE: FUN should return FUN(X) and not the sum-of-squares 
%   sum(FUN(X).^2)). (FUN(X) is summed and squared implicitly in the
%   algorithm.) 
%
%   X = LSQNONLIN(FUN,X0,LB,UB) defines a set of lower and upper bounds on
%   the design variables, X, so that the solution is in the range LB <= X
%   <= UB. Use empty matrices for LB and UB if no bounds exist. Set LB(i)
%   = -Inf if X(i) is unbounded below; set UB(i) = Inf if X(i) is
%   unbounded above.
%
%   X = LSQNONLIN(FUN,X0,LB,UB,OPTIONS) minimizes with the default
%   optimization parameters replaced by values in the structure OPTIONS, an
%   argument created with the OPTIMSET function. See OPTIMSET for details.
%   Used options are Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%   FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType,
%   LevenbergMarquardt, MaxFunEvals, MaxIter, DiffMinChange and
%   DiffMaxChange, LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG,
%   TypicalX, PlotFcns, and OutputFcn. Use the Jacobian option to specify
%   that FUN also returns a second output argument J that is the Jacobian
%   matrix at the point X. If FUN returns a vector F of m components when X
%   has length n, then J is an m-by-n matrix where J(i,j) is the partial
%   derivative of F(i) with respect to x(j). (Note that the Jacobian J is
%   the transpose of the gradient of F.)
%
%   X = LSQNONLIN(PROBLEM) solves the non-linear least squares problem 
%   defined in PROBLEM. PROBLEM is a structure with the function FUN in 
%   PROBLEM.objective, the start point in PROBLEM.x0, the lower bounds in 
%   PROBLEM.lb, the upper bounds in PROBLEM.ub, the options structure in 
%   PROBLEM.options, and solver name 'lsqnonlin' in PROBLEM.solver. Use 
%   this syntax to solve at the command line a problem exported from 
%   OPTIMTOOL. The structure PROBLEM must have all the fields. 
%
%   [X,RESNORM] = LSQNONLIN(FUN,X0,...) returns 
%   the value of the squared 2-norm of the residual at X: sum(FUN(X).^2). 
%
%   [X,RESNORM,RESIDUAL] = LSQNONLIN(FUN,X0,...) returns the value of the 
%   residual at the solution X: RESIDUAL = FUN(X).
%
%   [X,RESNORM,RESIDUAL,EXITFLAG] = LSQNONLIN(FUN,X0,...) returns an 
%   EXITFLAG that describes the exit condition of LSQNONLIN. Possible 
%   values of EXITFLAG and the corresponding exit conditions are 
%
%     1  LSQNONLIN converged to a solution X.
%     2  Change in X smaller than the specified tolerance.
%     3  Change in the residual smaller than the specified tolerance.
%     4  Magnitude search direction smaller than the specified tolerance.
%     0  Maximum number of function evaluations or of iterations reached.
%    -1  Algorithm terminated by the output function.
%    -2  Bounds are inconsistent.
%    -4  Line search cannot sufficiently decrease the residual along the 
%         current search direction.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQNONLIN(FUN,X0,...) returns a 
%   structure OUTPUT with the number of iterations taken in
%   OUTPUT.iterations, the number of function evaluations in
%   OUTPUT.funcCount, the algorithm used in OUTPUT.algorithm, the number
%   of CG iterations (if used) in OUTPUT.cgiterations, the first-order
%   optimality (if used) in OUTPUT.firstorderopt, and the exit message in
%   OUTPUT.message.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQNONLIN(FUN,X0,...) 
%   returns the set of Lagrangian multipliers, LAMBDA, at the solution: 
%   LAMBDA.lower for LB and LAMBDA.upper for UB.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = LSQNONLIN(FUN,
%   X0,...) returns the Jacobian of FUN at X.   
%
%   Examples
%     FUN can be specified using @:
%        x = lsqnonlin(@myfun,[2 3 4])
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x);
%
%   FUN can also be an anonymous function:
%
%       x = lsqnonlin(@(x) sin(3*x),[1 4])
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the non-linear 
%   least squares problem given in the function myfun, which is 
%   parameterized by its second argument c. Here myfun is an M-file 
%   function such as
%
%       function F = myfun(x,c)
%       F = [ 2*x(1) - exp(c*x(1))
%             -x(1) - exp(c*x(2))
%             x(1) - x(2) ];
%
%   To solve the least squares problem for a specific value of c, first 
%   assign the value to c. Then create a one-argument anonymous function 
%   that captures that value of c and calls myfun with two arguments. 
%   Finally, pass this anonymous function to LSQNONLIN:
%
%       c = -1; % define parameter first
%       x = lsqnonlin(@(x) myfun(x,c),[1;1])
%
%   See also OPTIMSET, LSQCURVEFIT, FSOLVE, @, INLINE.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/10 21:49:59 $

%   The default algorithm when OPTIONS.LargeScale = 'off' is the 
%   Levenberg-Marquardt method with a mixed quadratic and cubic line search
%   procedure.  A Gauss-Newton method is selected by setting 
%   OPTIONS.LargeScale = 'off' and OPTIONS.LevenbergMarquardt = 'off'. 
%
%

% ------------Initialization----------------


defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off','FunValCheck','off',...
   'Jacobian','off','JacobMult',[],... % JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',Inf,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','LevenbergMarquardt','on', ...
   'OutputFcn',[],'PlotFcns',[]); 

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 3, LB=[]; end
if nargin < 4, UB=[]; end
if nargin < 5, options=[]; end

if nargin == 1
    if isa(FUN,'struct')
        [FUN,x,LB,UB,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error('optim:lsqnonlin:InputArg','The input to LSQNONLIN should be either a structure with valid fields or consist of at least two arguments.');
    end
end

if nargin < 1 
  error('optim:lsqnonlin:NotEnoughInputs', ...
        'LSQNONLIN requires two input arguments.')
end

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(x,LB,UB);
    if ~isequal('double', dataType)
        error('optim:lsqnonlin:NonDoubleInput', ...
            'LSQNONLIN only accepts inputs of data type double.')
    end
catch
    error('optim:lsqnonlin:NonDoubleInput', ...
        'LSQNONLIN only accepts inputs of data type double.')
end

if nargout > 5
   computeLambda = 1;
else 
   computeLambda = 0;
end

caller = 'lsqnonlin'; XDATA = []; YDATA = [];
[x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = ...
   lsqncommon(FUN,x,XDATA,YDATA,LB,UB,options,defaultopt,caller,...
              computeLambda,varargin{:});
