function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%FMINCON finds a constrained minimum of a function of several variables.
%   FMINCON attempts to solve problems of the form:
%    min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%     X                     C(X) <= 0, Ceq(X) = 0   (nonlinear constraints)
%                              LB <= X <= UB        (bounds)
%                                                           
%   X = FMINCON(FUN,X0,A,B) starts at X0 and finds a minimum X to the 
%   function FUN, subject to the linear inequalities A*X <= B. FUN accepts 
%   input X and returns a scalar function value F evaluated at X. X0 may be
%   a scalar, vector, or matrix. 
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq) minimizes FUN subject to the linear 
%   equalities Aeq*X = Beq as well as A*X <= B. (Set A=[] and B=[] if no 
%   inequalities exist.)
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that a solution is found in 
%   the range LB <= X <= UB. Use empty matrices for LB and UB
%   if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; 
%   set UB(i) = Inf if X(i) is unbounded above.
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) subjects the minimization
%   to the constraints defined in NONLCON. The function NONLCON accepts X 
%   and returns the vectors C and Ceq, representing the nonlinear 
%   inequalities and equalities respectively. FMINCON minimizes FUN such 
%   that C(X) <= 0 and Ceq(X) = 0. (Set LB = [] and/or UB = [] if no bounds
%   exist.)
%
%   X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with 
%   the default optimization parameters replaced by values in the structure
%   OPTIONS, an argument created with the OPTIMSET function. See OPTIMSET
%   for details. For a list of options accepted by FMINCON refer to the
%   documentation.
%  
%   X = FMINCON(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the linear inequality constraints in PROBLEM.Aineq
%   and PROBLEM.bineq, the linear equality constraints in PROBLEM.Aeq and
%   PROBLEM.beq, the lower bounds in PROBLEM.lb, the upper bounds in 
%   PROBLEM.ub, the nonlinear constraint function in PROBLEM.nonlcon, the
%   options structure in PROBLEM.options, and solver name 'fmincon' in
%   PROBLEM.solver. Use this syntax to solve at the command line a problem 
%   exported from OPTIMTOOL. The structure PROBLEM must have all the fields.
%
%   [X,FVAL] = FMINCON(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG] = FMINCON(FUN,X0,...) returns an EXITFLAG that 
%   describes the exit condition of FMINCON. Possible values of EXITFLAG 
%   and the corresponding exit conditions are listed below.
%   
%   All algorithms:
%     1  First order optimality conditions satisfied to the specified 
%         tolerance.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Optimization terminated by the output function.
%    -2  No feasible point found.
%   Trust-region-reflective and interior-point:
%     2  Change in X less than the specified tolerance.
%   Trust-region-reflective:
%     3  Change in the objective function value less than the specified 
%         tolerance.
%   Active-set only:
%     4  Magnitude of search direction smaller than the specified tolerance
%         and constraint violation less than options.TolCon.
%     5  Magnitude of directional derivative less than the specified 
%         tolerance and constraint violation less than options.TolCon.
%   Interior-point:
%    -3  Problem appears to be unbounded.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINCON(FUN,X0,...) returns a structure 
%   OUTPUT with information such as total number of iterations, and final 
%   objective function value. See the documentation for a complete list.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINCON(FUN,X0,...) returns the 
%   Lagrange multipliers at the solution X: LAMBDA.lower for LB, 
%   LAMBDA.upper for UB, LAMBDA.ineqlin is for the linear inequalities, 
%   LAMBDA.eqlin is for the linear equalities, LAMBDA.ineqnonlin is for the
%   nonlinear inequalities, and LAMBDA.eqnonlin is for the nonlinear 
%   equalities.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD] = FMINCON(FUN,X0,...) returns the 
%   value of the gradient of FUN at the solution X.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = FMINCON(FUN,X0,...) 
%   returns the value of the exact or approximate Hessian of the Lagrangian
%   at X. 
%
%   Examples
%     FUN can be specified using @:
%        X = fmincon(@humps,...)
%     In this case, F = humps(X) returns the scalar function value F of 
%     the HUMPS function evaluated at X.
%
%     FUN can also be an anonymous function:
%        X = fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
%     returns X = [0;0].
%
%   If FUN or NONLCON are parameterized, you can use anonymous functions to
%   capture the problem-dependent parameters. Suppose you want to minimize 
%   the objective given in the function myfun, subject to the nonlinear 
%   constraint mycon, where these two functions are parameterized by their 
%   second argument a1 and a2, respectively. Here myfun and mycon are 
%   M-file functions such as
%
%        function f = myfun(x,a1)      
%        f = x(1)^2 + a1*x(2)^2;       
%                                      
%        function [c,ceq] = mycon(x,a2)
%        c = a2/x(1) - x(2);
%        ceq = [];
%
%   To optimize for specific values of a1 and a2, first assign the values 
%   to these two parameters. Then create two one-argument anonymous 
%   functions that capture the values of a1 and a2, and call myfun and 
%   mycon with two arguments. Finally, pass these anonymous functions to 
%   FMINCON:
%
%        a1 = 2; a2 = 1.5; % define parameters first
%        options = optimset('Algorithm','active-set'); % run active-set algorithm
%        x = fmincon(@(x) myfun(x,a1),[1;2],[],[],[],[],[],[],@(x) mycon(x,a2),options)
%
%   See also OPTIMSET, OPTIMTOOL, FMINUNC, FMINBND, FMINSEARCH, @, FUNCTION_HANDLE.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.7.4.1 $  $Date: 2008/01/21 15:53:49 $

defaultopt = struct( ...
    'Algorithm','trust-region-reflective', ...
    'AlwaysHonorConstraints','bounds', ...
    'DerivativeCheck','off', ...
    'Diagnostics','off', ...
    'DiffMaxChange',1e-1, ...
    'DiffMinChange',1e-8, ...
    'Display','final', ...
    'FinDiffType','forward', ...    
    'FunValCheck','off', ...
    'GradConstr','off', ...
    'GradObj','off', ...
    'HessFcn',[], ...
    'Hessian',[], ...    
    'HessMult',[], ...
    'HessPattern','sparse(ones(numberOfVariables))', ...
    'InitBarrierParam',0.1, ...
    'InitTrustRegionRadius','sqrt(numberOfVariables)', ...
    'LargeScale','on', ...
    'MaxFunEvals',[], ...
    'MaxIter',[], ...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
    'MaxProjCGIter','2*(numberOfVariables-numberOfEqualities)', ...    
    'MaxSQPIter','10*max(numberOfVariables,numberOfInequalities+numberOfBounds)', ...
    'NoStopIfFlatInfeas','off', ...
    'ObjectiveLimit',-1e20, ...
    'OutputFcn',[], ...
    'PhaseOneTotalScaling','off', ...
    'PlotFcns',[], ...
    'PrecondBandWidth',0, ...
    'RelLineSrchBnd',[], ...
    'RelLineSrchBndDuration',1, ...
    'ScaleProblem','obj-and-constr', ...
    'SubproblemAlgorithm','ldl-factorization', ...
    'TolCon',1e-6, ...
    'TolConSQP',1e-6, ...    
    'TolFun',1e-6, ...
    'TolGradCon',1e-6, ...
    'TolPCG',0.1, ...
    'TolProjCG',1e-2, ...
    'TolProjCGAbs',1e-10, ...
    'TolX',[], ...
    'TypicalX','ones(numberOfVariables,1)', ...
    'UseParallel','never' ...
    );

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
   X = defaultopt;
   return
end

if nargin < 10
    options = [];
    if nargin < 9
        NONLCON = [];
        if nargin < 8
            UB = [];
            if nargin < 7
                LB = [];
                if nargin < 6
                    Beq=[];
                    if nargin < 5
                        Aeq = [];
                    end
                end
            end
        end
    end
end

problemInput = false;
if nargin == 1
    if isa(FUN,'struct')
        problemInput = true;
        [FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error('optim:fmincon:InputArg','The input to FMINCON should be either a structure with valid fields or consist of at least four arguments.' );
    end
end

if nargin < 4 && ~problemInput
  error('optim:fmincon:AtLeastFourInputs','FMINCON requires at least four input arguments.')
end

if isempty(NONLCON) && isempty(A) && isempty(Aeq) && isempty(UB) && isempty(LB)
   error('optim:fmincon:ConstrainedProblemsOnly', ...
         'FMINCON is for constrained problems. Use FMINUNC for unconstrained problems.')
end

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(X,A,B,Aeq,Beq,LB,UB);
    if ~isequal('double', dataType)
        error('optim:fmincon:NonDoubleInput', ...
            'FMINCON only accepts inputs of data type double.')
    end
catch
    error('optim:fmincon:NonDoubleInput', ...
        'FMINCON only accepts inputs of data type double.')
end

if nargout > 4
   computeLambda = true;
else 
   computeLambda = false;
end

activeSet = 'medium-scale: SQP, Quasi-Newton, line-search';
trustRegionReflective = 'large-scale: trust-region reflective Newton';
interiorPoint = 'interior-point';

XOUT = X(:);
numberOfVariables=length(XOUT);
% Check for empty X
if numberOfVariables == 0
   error('optim:fmincon:EmptyX','You must provide a non-empty starting point.');
end

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
        verbosity = 4;
    otherwise
        verbosity = 2;
end

% Set to column vectors
B = B(:);
Beq = Beq(:);

LargeScaleFlag = strcmpi(optimget(options,'LargeScale',defaultopt,'fast'),'on'); 
Algorithm = optimget(options,'Algorithm',defaultopt,'fast'); 
% Option needed for processing initial guess
AlwaysHonorConstraints = optimget(options,'AlwaysHonorConstraints',defaultopt,'fast'); 

% Determine algorithm user chose via options. (We need this now
% to set OUTPUT.algorithm in case of early termination due to 
% inconsistent bounds.) 
% This algorithm choice may be modified later when we check the 
% problem type and supplied derivatives
algChoiceOptsConflict = false;
if strcmpi(Algorithm,'active-set')
    OUTPUT.algorithm = activeSet;
elseif strcmpi(Algorithm,'interior-point')
    OUTPUT.algorithm = interiorPoint;
else % Algorithm = trust-region-reflective
    if LargeScaleFlag
        OUTPUT.algorithm = trustRegionReflective;
    else
        % Conflicting options Algorithm='trust-region-reflective' and
        % LargeScale='off'. Choose active-set algorithm.
        algChoiceOptsConflict = true; % warn later, not in early termination
        OUTPUT.algorithm = activeSet;
    end
end

[XOUT,l,u,msg] = checkbounds(XOUT,LB,UB,numberOfVariables);
if ~isempty(msg)
   EXITFLAG = -2;
   [FVAL,LAMBDA,GRAD,HESSIAN] = deal([]);
   
   OUTPUT.iterations = 0;
   OUTPUT.funcCount = 0;
   OUTPUT.stepsize = [];
   if strcmpi(OUTPUT.algorithm,activeSet)
       OUTPUT.lssteplength = [];
   else % trust-region-reflective, interior-point
       OUTPUT.cgiterations = [];
   end
   if strcmpi(OUTPUT.algorithm,interiorPoint) || strcmpi(OUTPUT.algorithm,activeSet)
       OUTPUT.constrviolation = [];
   end
   OUTPUT.firstorderopt = [];
   OUTPUT.message = msg;
   
   X(:) = XOUT;
   if verbosity > 0
      disp(msg)
   end
   return
end
lFinite = l(~isinf(l));
uFinite = u(~isinf(u));

% Create structure of flags and initial values, initialize merit function
% type and the original shape of X.
flags.meritFunction = 0;
initVals.xOrigShape = X;

mtxmpy = optimget(options,'HessMult',defaultopt,'fast');
if isequal(mtxmpy,'hmult')
   warning('optim:fmincon:HessMultNameClash', ...
           ['Potential function name clash with a Toolbox helper function:\n',...
            ' Use a name besides ''hmult'' for your HessMult function to\n',...
            '  avoid errors or unexpected results.']);
end

diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
funValCheck = strcmpi(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
flags.grad = strcmpi(optimget(options,'GradObj',defaultopt,'fast'),'on');

% Notice that defaultopt.Hessian = [], so the variable "hessian" can be empty
hessian = optimget(options,'Hessian',defaultopt,'fast'); 
% If calling trust-region-reflective with an unavailable Hessian option value, 
% issue informative error message
if strcmpi(OUTPUT.algorithm,trustRegionReflective) && ...
        ~( isempty(hessian) || strcmpi(hessian,'on') || strcmpi(hessian,'user-supplied') || ...
           strcmpi(hessian,'off') || strcmpi(hessian,'fin-diff-grads')  )
    error('optim:fmincon:BadTRReflectHessianValue', ...
        ['Value of Hessian option unavailable in trust-region-reflective algorithm. Possible\n' ...
         ' values are ''user-supplied'' and ''fin-diff-grads''.'])
end

if ~iscell(hessian) && ( strcmpi(hessian,'user-supplied') || strcmpi(hessian,'on') )
    flags.hess = true;
else
    flags.hess = false;
end

if isempty(NONLCON)
   flags.constr = false;
else
   flags.constr = true;
end

% Process objective function
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
   % constrflag in optimfcnchk set to false because we're checking the objective, not constraint
   funfcn = optimfcnchk(FUN,'fmincon',length(varargin),funValCheck,flags.grad,flags.hess,false,Algorithm);
else
   error('optim:fmincon:InvalidFUN', ...
         ['FUN must be a function handle;\n', ...
          ' or, FUN may be a cell array that contains function handles.']);
end

% Process constraint function
if flags.constr % NONLCON is non-empty
   flags.gradconst = strcmpi(optimget(options,'GradConstr',defaultopt,'fast'),'on');
   % hessflag in optimfcnchk set to false because hessian is never returned by nonlinear constraint 
   % function
   %
   % constrflag in optimfcnchk set to true because we're checking the constraints
   confcn = optimfcnchk(NONLCON,'fmincon',length(varargin),funValCheck,flags.gradconst,false,true);
else
   flags.gradconst = false; 
   confcn{1} = '';
end

[rowAeq,colAeq] = size(Aeq);

% Look at problem type and supplied derivatives, and determine if need to
% switch algorithm
if strcmpi(OUTPUT.algorithm,activeSet)
    if algChoiceOptsConflict
        % Active-set algorithm chosen as a result of conflicting options
        warning('optim:fmincon:NLPAlgLargeScaleConflict', ...
            ['Options LargeScale = ''off'' and Algorithm = ''trust-region-reflective'' conflict.\n' ...
            'Ignoring Algorithm and running active-set method. To run trust-region-reflective, set\n' ...
            'LargeScale = ''on''. To run active-set without this warning, use Algorithm = ''active-set''.'])
    end
    if issparse(Aeq) || issparse(A)
        warning('optim:fmincon:ConvertingToFull', ...
            'Cannot use sparse matrices with active-set method: converting to full.')
    end
    if flags.hess % conflicting options
        flags.hess = false;
        warning('optim:fmincon:HessianIgnored', ...
            ['Active-set method is a Quasi-Newton method and does not use analytic Hessian.\n' ...
            'Hessian flag in options will be ignored (user-supplied Hessian will not be used).']);
        if isequal(funfcn{1},'fungradhess')
            funfcn{1}='fungrad';
        elseif  isequal(funfcn{1},'fun_then_grad_then_hess')
            funfcn{1}='fun_then_grad';
        end
    end
elseif strcmpi(OUTPUT.algorithm,trustRegionReflective)
    if isempty(NONLCON) && isempty(A) && isempty(Aeq) && flags.grad
        % if only l and u then call sfminbx
    elseif isempty(NONLCON) && isempty(A) && isempty(lFinite) && isempty(uFinite) && flags.grad ...
            && colAeq > rowAeq
        % if only Aeq beq and Aeq has more columns than rows, then call sfminle
    else
        warning('optim:fmincon:SwitchingToMediumScale', ...
            ['Trust-region-reflective method does not currently solve this type of problem,\n' ...
            ' using active-set (line search) instead.'])
        if isequal(funfcn{1},'fungradhess')
            funfcn{1}='fungrad';
            warning('optim:fmincon:HessianIgnored', ...
                ['Active-set method is a Quasi-Newton method and does not use\n' ...
                'analytic Hessian. Hessian flag in options will be ignored.'])
        elseif  isequal(funfcn{1},'fun_then_grad_then_hess')
            funfcn{1}='fun_then_grad';
            warning('optim:fmincon:HessianIgnored', ...
                ['Active-set method is a Quasi-Newton method and does not use\n' ...
                'analytic Hessian. Hessian flag in options will be ignored.'])
        end
        flags.hess = false;
        OUTPUT.algorithm = activeSet; % switch to active-set
    end
end

lenvlb = length(l);
lenvub = length(u);

if strcmpi(OUTPUT.algorithm,activeSet)
   %
   % Ensure starting point lies within bounds
   %
   i=1:lenvlb;
   lindex = XOUT(i)<l(i);
   if any(lindex),
      XOUT(lindex)=l(lindex)+1e-4; 
   end
   i=1:lenvub;
   uindex = XOUT(i)>u(i);
   if any(uindex)
      XOUT(uindex)=u(uindex);
   end
   X(:) = XOUT;
elseif strcmpi(OUTPUT.algorithm,trustRegionReflective)
   %
   % If components of initial x not within bounds, set those components  
   % of initial point to a "box-centered" point
   %
   arg = (u >= 1e10); arg2 = (l <= -1e10);
   u(arg) = inf;
   l(arg2) = -inf;
   xinitOutOfBounds_idx = XOUT < l | XOUT > u;
   if any(xinitOutOfBounds_idx)
       XOUT = startx(u,l,XOUT,xinitOutOfBounds_idx);
       X(:) = XOUT;
   end
else % interior-point
    % Variables: fixed, finite lower bounds, finite upper bounds
    xIndices = classifyBoundsOnVars(l,u,numberOfVariables);

    % If honor bounds mode, then check that initial point strictly satisfies the
    % simple inequality bounds on the variables and exactly satisfies fixed variable
    % bounds.
    if strcmpi(AlwaysHonorConstraints,'bounds')
        violatedFixedBnds_idx = XOUT(xIndices.fixed) ~= l(xIndices.fixed);
        violatedLowerBnds_idx = XOUT(xIndices.finiteLb) <= l(xIndices.finiteLb);
        violatedUpperBnds_idx = XOUT(xIndices.finiteUb) >= u(xIndices.finiteUb);
        if any(violatedLowerBnds_idx) || any(violatedUpperBnds_idx) || any(violatedFixedBnds_idx)
            if verbosity >= 4
                warning('optim:fmincon:XInitNotStrictInsideBnds', ...
                    ['Initial point not strictly inside bounds and AlwaysHonorConstraints=''bounds'';\n', ...
                     'shifting initial point inside bounds.'])
            end
            XOUT = shiftInitPtToInterior(numberOfVariables,XOUT,l,u,Inf);
            X(:) = XOUT;
        end
    end
end

% Evaluate function
initVals.g = zeros(numberOfVariables,1);
HESSIAN = [];

switch funfcn{1}
case 'fun'
   try
      initVals.f = funfcn{3}(X,varargin{:});
   catch
     error('optim:fmincon:ObjectiveError', ...
            ['FMINCON cannot continue because user supplied objective function' ...
             ' failed with the following error:\n%s'], lasterr)
   end
case 'fungrad'
   try
      [initVals.f,initVals.g(:)] = funfcn{3}(X,varargin{:});
   catch 
      error('optim:fmincon:ObjectiveError', ...
           ['FMINCON cannot continue because user supplied objective function' ...
            ' failed with the following error:\n%s'], lasterr)
   end
case 'fungradhess'
   try
      [initVals.f,initVals.g(:),HESSIAN] = funfcn{3}(X,varargin{:});
   catch
     error('optim:fmincon:ObjectiveError', ...
            ['FMINCON cannot continue because user supplied objective function' ...
             ' failed with the following error:\n%s'], lasterr)
   end
case 'fun_then_grad'
   try
      initVals.f = funfcn{3}(X,varargin{:});
   catch
     error('optim:fmincon:ObjectiveError', ...
            ['FMINCON cannot continue because user supplied objective function' ...
             ' failed with the following error:\n%s'], lasterr)
   end
   try
      initVals.g(:) = funfcn{4}(X,varargin{:});
   catch
      error('optim:fmincon:GradError', ...
            ['FMINCON cannot continue because user supplied objective gradient function' ...
             ' failed with the following error:\n%s'], lasterr)
   end
case 'fun_then_grad_then_hess'
   try
      initVals.f = funfcn{3}(X,varargin{:});
   catch
      error('optim:fmincon:ObjectiveError', ...
            ['FMINCON cannot continue because user supplied objective function' ...
             ' failed with the following error:\n%s'], lasterr)
   end
   try
      initVals.g(:) = funfcn{4}(X,varargin{:});
   catch
      error('optim:fmincon:GradientError', ...
            ['FMINCON cannot continue because user supplied objective gradient function' ...
             ' failed with the following error:\n%s'], lasterr)     
   end
   try
      HESSIAN = funfcn{5}(X,varargin{:});
   catch 
      error('optim:fmincon:HessianError', ...
            ['FMINCON cannot continue because user supplied objective Hessian function' ...
             ' failed with the following error:\n%s'], lasterr)     
   end
otherwise
   error('optim:fmincon:UndefinedCallType','Undefined calltype in FMINCON.');
end

% Check that the objective value is a scalar
if numel(initVals.f) ~= 1
   error('optim:fmincon:NonScalarObj','User supplied objective function must return a scalar value.')
end

% Evaluate constraints
switch confcn{1}
case 'fun'
   try 
      [ctmp,ceqtmp] = confcn{3}(X,varargin{:});
      initVals.ncineq = ctmp(:);
      initVals.nceq = ceqtmp(:);
      initVals.gnc = zeros(numberOfVariables,length(initVals.ncineq));
      initVals.gnceq = zeros(numberOfVariables,length(initVals.nceq));
   catch
      if findstr(xlate('Too many output arguments'),lasterr)
          if isa(confcn{3},'inline')
              error('optim:fmincon:InvalidInlineNonlcon', ...
                ['The inline function %s representing the constraints\n' ...
                 ' must return two outputs: the nonlinear inequality constraints and\n' ...
                 ' the nonlinear equality constraints.  At this time, inline objects may\n' ...
                 ' only return one output argument: use an M-file function instead.'], ...
                    formula(confcn{3}))            
          elseif isa(confcn{3},'function_handle')
              error('optim:fmincon:InvalidHandleNonlcon', ...
                   ['The constraint function %s must return two outputs:\n' ...
                    ' the nonlinear inequality constraints and\n' ...
                    ' the nonlinear equality constraints.'],func2str(confcn{3}))            
          else
              error('optim:fmincon:InvalidFunctionNonlcon', ...
                   ['The constraint function %s must return two outputs:\n' ...
                    ' the nonlinear inequality constraints and\n' ...
                    ' the nonlinear equality constraints.'],confcn{3})
          end
      else
        error('optim:fmincon:NonlconError', ... 
            ['FMINCON cannot continue because user supplied nonlinear constraint function\n' ...
            ' failed with the following error:\n%s'],lasterr)        
      end
   end
   
case 'fungrad'
   try
      [ctmp,ceqtmp,initVals.gnc,initVals.gnceq] = confcn{3}(X,varargin{:});
      initVals.ncineq = ctmp(:);
      initVals.nceq = ceqtmp(:);
   catch
      error('optim:fmincon:NonlconError', ... 
           ['FMINCON cannot continue because user supplied nonlinear constraint function\n' ...
            ' failed with the following error:\n%s'],lasterr)  
   end
case 'fun_then_grad'
   try
      [ctmp,ceqtmp] = confcn{3}(X,varargin{:});
      initVals.ncineq = ctmp(:);
      initVals.nceq = ceqtmp(:);
      [initVals.gnc,initVals.gnceq] = confcn{4}(X,varargin{:});
   catch
      error('optim:fmincon:NonlconFunOrGradError', ... 
           ['FMINCON cannot continue because user supplied nonlinear constraint function\n' ...
            'or nonlinear constraint gradient function failed with the following error:\n%s'],lasterr) 
   end
case ''
   initVals.ncineq=[];
   initVals.nceq =[];
   initVals.gnc = zeros(numberOfVariables,length(initVals.ncineq));
   initVals.gnceq = zeros(numberOfVariables,length(initVals.nceq));
otherwise
   error('optim:fmincon:UndefinedCalltype','Undefined calltype in FMINCON.');
end

non_eq = length(initVals.nceq);
non_ineq = length(initVals.ncineq);

% Make sure empty constraint and their derivatives have correct sizes (not 0-by-0):
if isempty(initVals.ncineq)
    initVals.ncineq = reshape(initVals.ncineq,0,1);
end
if isempty(initVals.nceq)
    initVals.nceq = reshape(initVals.nceq,0,1);
end
if isempty(Aeq)
    Aeq = reshape(Aeq,0,numberOfVariables);
    Beq = reshape(Beq,0,1);
end
if isempty(A)
    A = reshape(A,0,numberOfVariables);
    B = reshape(B,0,1);    
end
if isempty(initVals.gnc)
    initVals.gnc = reshape(initVals.gnc,numberOfVariables,0);
end
if isempty(initVals.gnceq)
    initVals.gnceq = reshape(initVals.gnceq,numberOfVariables,0);
end

[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);
[cgrow,cgcol] = size(initVals.gnc);
[ceqgrow,ceqgcol] = size(initVals.gnceq);

% These sizes checks assume that empty matrices have already been made the correct size
if Aeqcol ~= numberOfVariables
   error('optim:fmincon:WrongNumberOfColumnsInAeq','Aeq must have %i column(s).',numberOfVariables)
end
if Acol ~= numberOfVariables
   error('optim:fmincon:WrongNumberOfColumnsInA','A must have %i column(s).',numberOfVariables)
end
if cgrow ~= numberOfVariables || cgcol ~= non_ineq
   error('optim:fmincon:WrongSizeGradNonlinIneq', ...
         'Gradient of nonlinear inequality constraints must have size %i-by-%i.', ...
         numberOfVariables,non_ineq)
end
if ceqgrow ~= numberOfVariables || ceqgcol ~= non_eq
   error('optim:fmincon:WrongSizeGradNonlinEq', ...
         'Gradient of nonlinear equality constraints must have size %i-by-%i.', ...
         numberOfVariables,non_eq)
end

if diagnostics > 0
   % Do diagnostics on information so far
   diagnose('fmincon',OUTPUT,flags.grad,flags.hess,flags.constr,flags.gradconst,...
      ~LargeScaleFlag,options,defaultopt,XOUT,non_eq,...
      non_ineq,lin_eq,lin_ineq,l,u,funfcn,confcn,initVals.f,initVals.g,HESSIAN, ...
      initVals.ncineq,initVals.nceq,initVals.gnc,initVals.gnceq);
end

% call algorithm
if strcmpi(OUTPUT.algorithm,activeSet) % active-set
    defaultopt.MaxIter = 400; defaultopt.MaxFunEvals = '100*numberofvariables'; defaultopt.TolX = 1e-6;
    defaultopt.Hessian = 'off';
    [X,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN]=...
        nlconst(funfcn,X,l,u,full(A),B,full(Aeq),Beq,confcn,options,defaultopt, ...
        verbosity,flags,initVals,varargin{:});
elseif strcmpi(OUTPUT.algorithm,trustRegionReflective) % trust-region-reflective
   if (isequal(funfcn{1}, 'fun_then_grad_then_hess') || isequal(funfcn{1}, 'fungradhess'))
      Hstr = [];
   elseif (isequal(funfcn{1}, 'fun_then_grad') || isequal(funfcn{1}, 'fungrad'))
      n = length(XOUT); 
      Hstr = optimget(options,'HessPattern',defaultopt,'fast');
      if ischar(Hstr) 
         if isequal(lower(Hstr),'sparse(ones(numberofvariables))')
            Hstr = sparse(ones(n));
         else
            error('optim:fmincon:InvalidHessPattern', ...
                  'Option ''HessPattern'' must be a matrix if not the default.')
         end
      end
   end
   
   defaultopt.MaxIter = 400; defaultopt.MaxFunEvals = '100*numberofvariables'; defaultopt.TolX = 1e-6;
   defaultopt.Hessian = 'off';
   if isempty(Aeq)
      [X,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
         sfminbx(funfcn,X,l,u,verbosity,options,defaultopt,computeLambda,initVals.f,initVals.g, ...
         HESSIAN,Hstr,varargin{:});
   else
      [X,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
         sfminle(funfcn,X,sparse(Aeq),Beq,verbosity,options,defaultopt,computeLambda,initVals.f, ...
         initVals.g,HESSIAN,Hstr,varargin{:});
   end
else % interior-point
    defaultopt.MaxIter = 1000; defaultopt.MaxFunEvals = 3000; defaultopt.TolX = 1e-10;
    defaultopt.Hessian = 'bfgs';
    mEq = lin_eq + non_eq + nnz(xIndices.fixed); % number of equalities
    % Interior-point-specific options. Default values for lbfgs memory is 10, and 
    % ldl pivot threshold is 0.01
    options_ip = getIpOptions(options,numberOfVariables,mEq,flags.constr,defaultopt,10,0.01); 
    [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = barrier(funfcn,X,A,B,Aeq,Beq,l,u,confcn,options_ip.HessFcn, ...
        initVals.f,initVals.g,initVals.ncineq,initVals.nceq,initVals.gnc,initVals.gnceq,HESSIAN, ...
        xIndices,options_ip,varargin{:});
end





