function [x,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(FUN,x,options,varargin)
%FMINUNC finds a local minimum of a function of several variables.
%   X = FMINUNC(FUN,X0) starts at X0 and attempts to find a local minimizer
%   X of the function FUN. FUN accepts input X and returns a scalar
%   function value F evaluated at X. X0 can be a scalar, vector or matrix. 
%
%   X = FMINUNC(FUN,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%   FunValCheck, GradObj, HessPattern, Hessian, HessMult, HessUpdate,
%   InitialHessType, InitialHessMatrix, MaxFunEvals, MaxIter, DiffMinChange
%   and DiffMaxChange, LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG,
%   PlotFcns, OutputFcn, and TypicalX. Use the GradObj option to specify
%   that FUN also returns a second output argument G that is the partial
%   derivatives of the function df/dX, at the point X. Use the Hessian
%   option to specify that FUN also returns a third output argument H that
%   is the 2nd partial derivatives of the function (the Hessian) at the
%   point X. The Hessian is only used by the large-scale method, not the
%   line-search method. 
%
%   X = FMINUNC(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminunc' in PROBLEM.solver. Use this syntax to solve at the 
%   command line a problem exported from OPTIMTOOL. The structure PROBLEM 
%   must have all the fields.
%
%   [X,FVAL] = FMINUNC(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = FMINUNC(FUN,X0,...) returns an EXITFLAG that 
%   describes the exit condition of FMINUNC. Possible values of EXITFLAG 
%   and the corresponding exit conditions are
%
%     1  Magnitude of gradient smaller than the specified tolerance. 
%     2  Change in X smaller than the specified tolerance.
%     3  Change in the objective function value smaller than the specified 
%         tolerance (only occurs in the large-scale method).
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Algorithm terminated by the output function.
%    -2  Line search cannot find an acceptable point along the current
%         search direction (only occurs in the medium-scale method).
%   
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINUNC(FUN,X0,...) returns a structure 
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the 
%   number of function evaluations in OUTPUT.funcCount, the algorithm used 
%   in OUTPUT.algorithm, the number of CG iterations (if used) in
%   OUTPUT.cgiterations, the first-order optimality (if used) in
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINUNC(FUN,X0,...) returns the value 
%   of the gradient of FUN at the solution X.
%
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = FMINUNC(FUN,X0,...) returns the 
%   value of the Hessian of the objective function FUN at the solution X.
%
%   Examples
%     FUN can be specified using @:
%        X = fminunc(@myfun,2)
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x) + 3;
%
%     To minimize this function with the gradient provided, modify
%     the function myfun so the gradient is the second output argument:
%        function [f,g] = myfun(x)
%         f = sin(x) + 3;
%         g = cos(x);
%     and indicate the gradient value is available by creating an options
%     structure with OPTIONS.GradObj set to 'on' (using OPTIMSET):
%        options = optimset('GradObj','on');
%        x = fminunc(@myfun,4,options);
%
%     FUN can also be an anonymous function:
%        x = fminunc(@(x) 5*x(1)^2 + x(2)^2,[5;1])
%
%   If FUN is parameterized, you can use anonymous functions to capture the
%   problem-dependent parameters. Suppose you want to minimize the 
%   objective given in the function myfun, which is parameterized by its 
%   second argument c. Here myfun is an M-file function such as
%
%     function [f,g] = myfun(x,c)
%
%     f = c*x(1)^2 + 2*x(1)*x(2) + x(2)^2; % function
%     g = [2*c*x(1) + 2*x(2)               % gradient
%          2*x(1) + 2*x(2)];
%
%   To optimize for a specific value of c, first assign the value to c. 
%   Then create a one-argument anonymous function that captures that value 
%   of c and calls myfun with two arguments. Finally, pass this anonymous 
%   function to FMINUNC:
%
%     c = 3;                              % define parameter first
%     options = optimset('GradObj','on'); % indicate gradient is provided 
%     x = fminunc(@(x) myfun(x,c),[1;1],options)
%
%   See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, @, INLINE.

%   When options.LargeScale=='on', the algorithm is a trust-region method.
%   When options.LargeScale=='off', the algorithm is the BFGS Quasi-Newton 
%   method with a mixed quadratic and cubic line search procedure. 

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/10 21:49:53 $
 
% ------------Initialization----------------
defaultopt = struct('Display','final','LargeScale','on', ...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...   
   'Diagnostics','off','FunValCheck','off',...
   'GradObj','off','MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'Hessian','off','HessMult',[],...
   'HessPattern','sparse(ones(numberOfVariables))',...
   'HessUpdate','bfgs','OutputFcn',[],'PlotFcns',[], ...
   'InitialHessType','scaled-identity','InitialHessMatrix',[]); 

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 3, options=[]; end 

% Detect problem structure input
if nargin == 1
    if isa(FUN,'struct')
        [FUN,x,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error('optim:fminunc:InputArg','The input to FMINUNC should be either a structure with valid fields or consist of at least two arguments.');
    end
end

if nargin == 0 
  error('optim:fminunc:NotEnoughInputs','FMINUNC requires two input arguments.')
end

if nargout > 5
  computeHessian = true;
else
  computeHessian = false;    
end

% Check for non-double inputs
if ~isa(x,'double')
  error('optim:fminunc:NonDoubleInput', ...
        'FMINUNC only accepts inputs of data type double.')
end

XOUT=x(:);
numberOfVariables=length(XOUT);
medium = 'medium-scale: Quasi-Newton line search'; 
large = 'large-scale: trust-region Newton'; 

switch optimget(options,'Display',defaultopt,'fast')
case {'off','none'}
   verbosity = 0;
case 'notify'
   verbosity = 1;  
case 'final'
   verbosity = 2;   
case 'iter'
   verbosity = 3;
case 'testing'
   verbosity = Inf;
otherwise
   verbosity = 2;
end
diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
mtxmpy = optimget(options,'HessMult',defaultopt,'fast');
if isequal(mtxmpy,'hmult')
   warning('optim:fminunc:HessMultNameClash', ...
      ['Potential function name clash with a Toolbox helper function:\n' ...
       ' Use a name besides ''hmult'' for your HessMult function to' ...
       ' avoid errors\n or unexpected results.'])
end
gradflag =  strcmp(optimget(options,'GradObj',defaultopt,'fast'),'on');
Hessian = optimget(options,'Hessian',defaultopt,'fast');

% line_search: 0 means trust-region, 1 means line-search
line_search = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off');

if ( strcmpi(Hessian,'on') || strcmpi(Hessian,'user-supplied') )
    hessflag = true;
elseif strcmpi(Hessian,'off') || strcmpi(Hessian,'fin-diff-grads')
    hessflag = false;
else
    % If calling large-scale algorithm with an unavailable Hessian option value,
    % issue informative error message
    if ~line_search
        error('optim:fminunc:BadTRReflectHessianValue', ...
            ['Value of Hessian option unavailable in large-scale algorithm. Possible\n' ...
            ' values are ''on'' and ''off''.'])
    end
end

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
computeLambda = 0;

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
   funfcn = optimfcnchk(FUN,'fminunc',length(varargin),funValCheck,gradflag,hessflag);
else
   error('optim:fminunc:InvalidFUN','FUN must be function handle or a cell array of two function handles.')
end

GRAD = zeros(numberOfVariables,1);
HESS = [];

switch funfcn{1}
case 'fun'
   f = feval(funfcn{3},x,varargin{:});
case 'fungrad'
   [f,GRAD(:)] = feval(funfcn{3},x,varargin{:});
case 'fungradhess'
   [f,GRAD(:),HESS] = feval(funfcn{3},x,varargin{:});
case 'fun_then_grad'
   f = feval(funfcn{3},x,varargin{:}); 
   GRAD(:) = feval(funfcn{4},x,varargin{:});
case 'fun_then_grad_then_hess'
   f = feval(funfcn{3},x,varargin{:}); 
   GRAD(:) = feval(funfcn{4},x,varargin{:});
   HESS = feval(funfcn{5},x,varargin{:});
otherwise
   error('optim:fminunc:UndefCalltype','Undefined calltype in FMINUNC.');
end

% Check that the objective value is a scalar
if numel(f) ~= 1
   error('optim:fminunc:NonScalarObj','User supplied objective function must return a scalar value.')
end

% Determine algorithm
% If line-search and no hessian,  then call line-search algorithm
if line_search  && ...
      (~isequal(funfcn{1}, 'fun_then_grad_then_hess') && ~isequal(funfcn{1}, 'fungradhess'))
  output.algorithm = medium; 
    
  % Line-search and Hessian -- no can do, so do line-search after warning: ignoring hessian.   
elseif line_search && ...
        (isequal(funfcn{1}, 'fun_then_grad_then_hess') || isequal(funfcn{1}, 'fungradhess'))
    warning('optim:fminunc:HessIgnored', ...  
        ['Medium-scale method is a Quasi-Newton method and does not use analytic Hessian.\n' ...
        ' Hessian flag in options will be ignored (user-supplied Hessian will not be used).'])
    if isequal(funfcn{1}, 'fun_then_grad_then_hess')
        funfcn{1} = 'fun_then_grad';
    elseif isequal(funfcn{1}, 'fungradhess')
        funfcn{1} = 'fungrad';
    end
    output.algorithm = medium;
    % If not line-search (trust-region) and Hessian, call trust-region   
elseif ~line_search && ...
        (isequal(funfcn{1}, 'fun_then_grad_then_hess') || isequal(funfcn{1}, 'fungradhess'))
   l=[]; u=[]; Hstr=[];
   output.algorithm = large; 
% If not line search (trust-region) and no Hessian but grad, use sparse finite-differencing.
elseif ~line_search && ...
      (isequal(funfcn{1}, 'fun_then_grad') || isequal(funfcn{1}, 'fungrad'))
   n = length(XOUT); 
   Hstr = optimget(options,'HessPattern',defaultopt,'fast');
   if ischar(Hstr) 
      if isequal(lower(Hstr),'sparse(ones(numberofvariables))')
      % Put this code separate as it might generate OUT OF MEMORY error
         Hstr = sparse(ones(n));
      else
         error('optim:fminunc:InvalidHessPattern', ...
               'Option ''HessPattern'' must be a matrix if not the default.')
      end
   end
   l=[]; u=[];
   output.algorithm = large;
   
   % Trust region but no grad, no can do; warn and use line-search    
elseif ~line_search
   warning('optim:fminunc:SwitchingMethod', ...
      ['Gradient must be provided for trust-region method;\n ' ...
      ' using line-search method instead.'])
   output.algorithm = medium;
else
   error('optim:fminunc:InvalidProblem','Problem not handled by FMINUNC.')   
end

if diagnostics > 0
   % Do diagnostics on information so far
   constflag = 0; gradconstflag = 0; non_eq=0;non_ineq=0;lin_eq=0;lin_ineq=0;
   LB=[]; UB =[];confcn{1}=[];c=[];ceq=[];cGRAD=[];ceqGRAD=[];
   diagnose('fminunc',output,gradflag,hessflag,constflag,gradconstflag,...
      line_search,options,defaultopt,XOUT,non_eq,...
      non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn,f,GRAD,HESS,c,ceq,cGRAD,ceqGRAD);
   
end

% If line-search and no hessian,  then call line-search algorithm
if isequal(output.algorithm, medium)
   [x,FVAL,GRAD,HESSIAN,EXITFLAG,OUTPUT] = fminusub(funfcn,x,verbosity, ...
      options,defaultopt,f,GRAD,HESS,computeHessian,varargin{:});
elseif isequal(output.algorithm, large)
   [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = sfminbx(funfcn,x,l,u, ...
      verbosity,options,defaultopt,computeLambda,f,GRAD,HESS,Hstr,varargin{:});
   OUTPUT.algorithm = large; % override sfminbx output: not using the reflective 
                             % part of the method   
end


