function [x,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(FUN,x,options,varargin)
%FSOLVE solves systems of nonlinear equations of several variables.
%
%   FSOLVE attempts to solve equations of the form:
%             
%   F(X) = 0    where F and X may be vectors or matrices.   
%
%   X = FSOLVE(FUN,X0) starts at the matrix X0 and tries to solve the 
%   equations in FUN.  FUN accepts input X and returns a vector (matrix) of 
%   equation values F evaluated at X. 
%
%   X = FSOLVE(FUN,X0,OPTIONS) solves the equations with the default 
%   optimization parameters replaced by values in the structure OPTIONS, an
%   argument created with the OPTIMSET function.  See OPTIMSET for details.
%   Used options are Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%   FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType,
%   NonlEqnAlgorithm, MaxFunEvals, MaxIter, PlotFcns, OutputFcn,
%   DiffMinChange and DiffMaxChange, LargeScale, MaxPCGIter,
%   PrecondBandWidth, TolPCG, and TypicalX. Use the Jacobian option to
%   specify that FUN also returns a second output argument J that is the
%   Jacobian matrix at the point X. If FUN returns a vector F of m
%   components when X has length n, then J is an m-by-n matrix where J(i,j)
%   is the partial derivative of F(i) with respect to x(j). (Note that the
%   Jacobian J is the transpose of the gradient of F.)
%
%   X = FSOLVE(PROBLEM) solves system defined in PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fsolve' in PROBLEM.solver.  Use this syntax to solve at the 
%   command line a problem exported from OPTIMTOOL. The structure PROBLEM 
%   must have all the fields.
%
%   [X,FVAL] = FSOLVE(FUN,X0,...) returns the value of the equations FUN 
%   at X. 
%
%   [X,FVAL,EXITFLAG] = FSOLVE(FUN,X0,...) returns an EXITFLAG that 
%   describes the exit condition of FSOLVE. Possible values of EXITFLAG and
%   the corresponding exit conditions are
%
%     1  FSOLVE converged to a solution X.
%     2  Change in X smaller than the specified tolerance.
%     3  Change in the residual smaller than the specified tolerance.
%     4  Magnitude of search direction smaller than the specified 
%         tolerance.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Algorithm terminated by the output function.
%    -2  Algorithm seems to be converging to a point that is not a root.
%    -3  Trust region radius became too small.
%    -4  Line search cannot sufficiently decrease the residual along the 
%         current search direction.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FSOLVE(FUN,X0,...) returns a structure 
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the 
%   number of function evaluations in OUTPUT.funcCount, the algorithm used 
%   in OUTPUT.algorithm, the number of CG iterations (if used) in 
%   OUTPUT.cgiterations, the first-order optimality (if used) in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,JACOB] = FSOLVE(FUN,X0,...) returns the 
%   Jacobian of FUN at X.  
%
%   Examples
%     FUN can be specified using @:
%        x = fsolve(@myfun,[2 3 4],optimset('Display','iter'))
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x);
%
%   FUN can also be an anonymous function:
%
%       x = fsolve(@(x) sin(3*x),[1 4],optimset('Display','off'))
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the system of 
%   nonlinear equations given in the function myfun, which is parameterized 
%   by its second argument c. Here myfun is an M-file function such as
%     
%       function F = myfun(x,c)
%       F = [ 2*x(1) - x(2) - exp(c*x(1))
%             -x(1) + 2*x(2) - exp(c*x(2))];
%           
%   To solve the system of equations for a specific value of c, first 
%   assign the value to c. Then create a one-argument anonymous function 
%   that captures that value of c and calls myfun with two arguments. 
%   Finally, pass this anonymous function to FSOLVE:
%
%       c = -1; % define parameter first
%       x = fsolve(@(x) myfun(x,c),[-5;-5])
%
%   See also OPTIMSET, LSQNONLIN, @, INLINE.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/12/10 21:49:55 $

% ------------Initialization----------------

defaultopt = struct('Display','final','LargeScale','off',...
   'NonlEqnAlgorithm','dogleg',...
   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
   'Diagnostics','off','FunValCheck','off',...
   'Jacobian','off','JacobMult',[],...% JacobMult set to [] by default
   'JacobPattern','sparse(ones(Jrows,Jcols))',...
   'MaxFunEvals','100*numberOfVariables',...
   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
   'PrecondBandWidth',Inf,'TypicalX','ones(numberOfVariables,1)',...
   'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
   'TolPCG',0.1,'MaxIter',400,...
   'LineSearchType','quadcubic','OutputFcn',[],'PlotFcns',[]);

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
        error('optim:fsolve:InputArg','The input to FSOLVE should be either a structure with valid fields or consist of at least two arguments.');
    end
end

if nargin == 0
  error('optim:fsolve:NotEnoughInputs','FSOLVE requires at least two input arguments.')
end

% Check for non-double inputs
if ~isa(x,'double')
  error('optim:fsolve:NonDoubleInput', ...
        'FSOLVE only accepts inputs of data type double.')
end

LB = []; UB = []; 
xstart=x(:);
numberOfVariables=length(xstart);

large        = 'large-scale';
medium       = 'medium-scale: line search';
dogleg       = 'trust-region dogleg';

switch optimget(options,'Display',defaultopt,'fast')
    case {'off','none'}
        verbosity = 0;
    case 'iter'
        verbosity = 2;
    case 'final'
        verbosity = 1;
    case 'testing'
        verbosity = Inf;
    otherwise
        verbosity = 1;
end
diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
% 0 means large-scale trust-region, 1 means medium-scale algorithm
mediumflag = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
switch optimget(options,'NonlEqnAlgorithm',defaultopt,'fast')
    case 'dogleg'
        algorithmflag = 1;
    case 'lm'
        algorithmflag = 2;
    case 'gn'
        algorithmflag = 3;
    otherwise
        algorithmflag = 1;
end
mtxmpy = optimget(options,'JacobMult',defaultopt,'fast');
if isequal(mtxmpy,'atamult')
    warning('optim:fsolve:NameClash', ...
        ['Potential function name clash with a Toolbox helper function:\n' ...
        'Use a name besides ''atamult'' for your JacobMult function to\n' ...
        'avoid errors or unexpected results.'])
end

% Convert to inline function as needed
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    funfcn = lsqfcnchk(FUN,'fsolve',length(varargin),funValCheck,gradflag);
else
    error('optim:fsolve:InvalidFUN', ...
        ['FUN must be a function name, valid string expression, or inline object;\n' ...
        ' or, FUN may be a cell array that contains these type of objects.'])
end

JAC = [];
x(:) = xstart;
switch funfcn{1}
    case 'fun'
        fuser = feval(funfcn{3},x,varargin{:});
        f = fuser(:);
        nfun=length(f);
    case 'fungrad'
        [fuser,JAC] = feval(funfcn{3},x,varargin{:});
        f = fuser(:);
        nfun=length(f);
    case 'fun_then_grad'
        fuser = feval(funfcn{3},x,varargin{:});
        f = fuser(:);
        JAC = feval(funfcn{4},x,varargin{:});
        nfun=length(f);
    otherwise
        error('optim:fsolve:UndefinedCalltype','Undefined calltype in FSOLVE.')
end

if gradflag
    % check size of JAC
    [Jrows, Jcols]=size(JAC);
    if isempty(mtxmpy)
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows~=nfun || Jcols~=numberOfVariables
            error('optim:fsolve:InvalidJacobian', ...
                ['User-defined Jacobian is not the correct size:\n' ...
                ' the Jacobian matrix should be %d-by-%d.'],nfun,numberOfVariables)
        end
    end
else
    Jrows = nfun;
    Jcols = numberOfVariables;
end

XDATA = []; YDATA = []; caller = 'fsolve';

% Choose what algorithm to run: determine (i) OUTPUT.algorithm and 
% (ii) if and only if OUTPUT.algorithm = medium, also option.LevenbergMarquardt.
% Option LevenbergMarquardt is used internally; it's not user settable. For
% this reason we change this option directly, for speed; users should use
% optimset.
if ~mediumflag 
    if nfun >= numberOfVariables
        % large-scale method and enough equations (as many as variables)
        OUTPUT.algorithm = large;
    else 
        % large-scale method and not enough equations - switch to medium-scale algorithm
        warning('optim:fsolve:FewerFunsThanVars', ...
                ['Large-scale method requires at least as many equations as variables;\n' ...
                 ' using line-search method instead.'])
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'off'; 
    end
else
    if algorithmflag == 1 && nfun == numberOfVariables 
        OUTPUT.algorithm = dogleg;
    elseif algorithmflag == 1 && nfun ~= numberOfVariables
        warning('optim:fsolve:NonSquareSystem', ...
                ['Default trust-region dogleg method of FSOLVE cannot\n handle non-square systems; ', ...
                 'using Gauss-Newton method instead.']);
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'off';
    elseif algorithmflag == 2
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'on';
    else % algorithmflag == 3
        OUTPUT.algorithm = medium;
        options.LevenbergMarquardt = 'off';
    end
end

if diagnostics > 0
    % Do diagnostics on information so far
    constflag = 0; gradconstflag = 0; non_eq=0;non_ineq=0;lin_eq=0;lin_ineq=0;
    confcn{1}=[];c=[];ceq=[];cGRAD=[];ceqGRAD=[];
    hessflag = 0; HESS=[];
    diagnose('fsolve',OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
        mediumflag,options,defaultopt,xstart,non_eq,...
        non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn,f,JAC,HESS,c,ceq,cGRAD,ceqGRAD);

end

% Execute algorithm
if isequal(OUTPUT.algorithm, large)
    if ~gradflag
        Jstr = optimget(options,'JacobPattern',defaultopt,'fast');
        if ischar(Jstr)
            if isequal(lower(Jstr),'sparse(ones(jrows,jcols))')
                Jstr = sparse(ones(Jrows,Jcols));
            else
                error('optim:fsolve:InvalidJacobPattern', ...
                    'Option ''JacobPattern'' must be a matrix if not the default.')
            end
        end
    else
        Jstr = [];
    end
    computeLambda = 0;
    [x,FVAL,LAMBDA,JACOB,EXITFLAG,OUTPUT,msg]=...
        snls(funfcn,x,LB,UB,verbosity,options,defaultopt,f,JAC,XDATA,YDATA,caller,...
        Jstr,computeLambda,varargin{:});
elseif isequal(OUTPUT.algorithm, dogleg)
    % trust region dogleg method
    Jstr = [];
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg]=...
        trustnleqn(funfcn,x,verbosity,gradflag,options,defaultopt,f,JAC,...
        Jstr,varargin{:});
else
    % line search (Gauss-Newton or Levenberg-Marquardt)
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msg] = ...
        nlsq(funfcn,x,verbosity,options,defaultopt,f,JAC,XDATA,YDATA,caller,varargin{:});
end

Resnorm = FVAL'*FVAL;  % assumes FVAL still a vector
if EXITFLAG > 0 % if we think we converged:
    if Resnorm > sqrt(optimget(options,'TolFun',defaultopt,'fast'))
        OUTPUT.message = ...
            sprintf(['Optimizer appears to be converging to a minimum that is not a root:\n' ...
            'Sum of squares of the function values is > sqrt(options.TolFun).\n' ...
            'Try again with a new starting point.']);
        if verbosity > 0
            disp(OUTPUT.message)
        end
        EXITFLAG = -2;
    else
        OUTPUT.message = msg;
        if verbosity > 0
            disp(OUTPUT.message);
        end
    end
else
    OUTPUT.message = msg;
    if verbosity > 0
        disp(OUTPUT.message);
    end
end

% Reset FVAL to shape of the user-function output, fuser
FVAL = reshape(FVAL,size(fuser));

