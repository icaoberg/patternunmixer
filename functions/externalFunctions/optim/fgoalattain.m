function [x,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT,LAMBDA] = fgoalattain(FUN,x,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%FGOALATTAIN solves the multi-objective goal attainment optimization 
% problem.
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT)
%   tries to make the objective functions (F) supplied by the function FUN
%   attain the goals (GOAL) by varying X. The goals are weighted according 
%   to WEIGHT. In doing so the following nonlinear programming problem is 
%   solved:
%            min     { LAMBDA :  F(X)-WEIGHT.*LAMBDA<=GOAL } 
%          X,LAMBDA  
%
%   FUN accepts input X and returns a vector (matrix) of function values F 
%   evaluated at X. X0 may be a scalar, vector, or matrix.  
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B) solves the goal attainment 
%   problem subject to the linear inequalities A*X <= B.
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq) solves the goal
%   attainment problem subject to the linear equalities Aeq*X = Beq as
%   well.  
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB) defines a set of 
%   lower and upper bounds on the design variables, X, so that the solution
%   is in the range LB <= X <= UB. Use empty matrices for LB and U if no 
%   bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; set 
%   UB(i) = Inf if X(i) is unbounded above.
%   
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON) subjects 
%   the goal attainment problem to the constraints defined in NONLCON 
%   (usually an M-file: NONLCON.m). The function NONLCON should return the 
%   vectors C and Ceq, representing the nonlinear inequalities and 
%   equalities respectively, when called with feval: 
%   [C, Ceq] = feval(NONLCON,X). FGOALATTAIN optimizes such that C(X) <= 0 
%   and Ceq(X) = 0.
%
%   X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
%   minimizes the with default optimization parameters replaced by values
%   in the structure OPTIONS, an argument created with the OPTIMSET
%   function.  See OPTIMSET for details. Used options are Display, TolX,
%   TolFun, TolCon, DerivativeCheck, FunValCheck, GradObj, GradConstr,
%   MaxFunEvals, MaxIter, MeritFunction, GoalsExactAchieve, Diagnostics,
%   DiffMinChange, DiffMaxChange, PlotFcns, OutputFcn, and TypicalX. Use
%   the GradObj option to specify that FUN may be called with two output
%   arguments where the second, G, is the partial derivatives of the
%   function df/dX, at the point X: [F,G] = feval(FUN,X). Use the
%   GradConstr option to specify that NONLCON may be called with four
%   output arguments: [C,Ceq,GC,GCeq] = feval(NONLCON,X) where GC is the
%   partial derivatives of the constraint vector of inequalities C an GCeq
%   is the partial derivatives of the constraint vector of equalities Ceq.
%   Use OPTIONS = [] as a place holder if no options are set.
%
%   X = FGOALATTAIN(PROBLEM) solves the goal attainment problem defined in 
%   PROBLEM. PROBLEM is a structure with the function FUN in 
%   PROBLEM.objective, the start point in PROBLEM.x0, the 'goal' vector in 
%   PROBLEM.goal, the 'weight' vector in PROBLEM.weight, the linear 
%   inequality constraints in PROBLEM.Aineq and PROBLEM.bineq, the linear 
%   equality constraints in PROBLEM.Aeq and PROBLEM.beq, the lower bounds 
%   in PROBLEM.lb, the upper bounds in PROBLEM.ub, the nonlinear constraint
%   function in PROBLEM.nonlcon, the options structure in PROBLEM.options, 
%   and solver name 'fgoalattain' in PROBLEM.solver. Use this syntax to 
%   solve at the command line a problem exported from OPTIMTOOL. The 
%   structure PROBLEM must have all the fields.
%
%   [X,FVAL] = FGOALATTAIN(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,ATTAINFACTOR] = FGOALATTAIN(FUN,X0,...) returns the attainment
%   factor at the solution X. If ATTAINFACTOR is negative, the goals have
%   been over- achieved; if ATTAINFACTOR is positive, the goals have been
%   under-achieved.
%
%   [X,FVAL,ATTAINFACTOR,EXITFLAG] = FGOALATTAIN(FUN,X0,...) returns an 
%   EXITFLAG that describes the exit condition of FGOALATTAIN. Possible 
%   values of EXITFLAG and the corresponding exit conditions are
%
%     1  FGOALATTAIN converged to a solution X.
%     4  Magnitude of search direction smaller than the specified tolerance
%         and constraint violation less than options.TolCon.
%     5  Magnitude of directional derivative less than the specified 
%         tolerance and constraint violation less than options.TolCon.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Optimization terminated by the output function.
%    -2  No feasible point found.
%   
%   [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT] = FGOALATTAIN(FUN,X0,...) returns 
%   a structure OUTPUT with the number of iterations taken in 
%   OUTPUT.iterations, the number of function evaluations in 
%   OUTPUT.funcCount, the norm of the final step in OUTPUT.stepsize, the 
%   final line search steplength in OUTPUT.lssteplength, the algorithm used
%   in OUTPUT.algorithm, the first-order optimality in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
% 
%   [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT,LAMBDA] = FGOALATTAIN(FUN,X0,...)
%   returns the Lagrange multiplier at the solution X: LAMBDA.lower for
%   LB, LAMBDA.upper for UB, LAMBDA.ineqlin is for the linear
%   inequalities, LAMBDA.eqlin is for the linear equalities,
%   LAMBDA.ineqnonlin is for the nonlinear inequalities, and
%   LAMBDA.eqnonlin is for the nonlinear equalities.
%
%   For more details, type the M-file FGOALATTAIN.M.
%
%   See also OPTIMSET, OPTIMGET.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.7.4.1 $  $Date: 2008/01/21 15:53:48 $

% ---------------------More Details---------------------------
% [x]=fgoalattain(F,x,GOAL,WEIGHT,[],[],[],[],[],[],[],OPTIONS)
% Solves the goal attainment problem where:
%
%  X  Is a set of design parameters which can be varied.
%  F  Is a set of objectives which are dependent on X.
%  GOAL Set of design goals. The optimizer will try to make 
%         F<GOAL, F=GOAL, or F>GOAL depending on the formulation.
%  WEIGHT Set of weighting parameters which determine the 
%         relative under or over achievement of the objectives.
%         Notes:
%           1.Setting WEIGHT=abs(GOAL)  will try to make the objectives
%             less than the goals resulting in roughly the same 
%             percentage under or over achievement of the goals.
%             Note: use WEIGHT 1 for GOALS that are 0 (see Note 3 below).
%           2. Setting WEIGHT=-abs(GOAL) will try to make the objectives
%              greater then the goals resulting in roughly the same percentage 
%              under- or over-achievement in the goals.
%             Note: use WEIGHT 1 for GOALS that are 0 (see Note 3 below).
%           3. Setting WEIGHT(i)=0  indicates a hard constraint.
%              i.e. F<=GOAL.
%  OPTIONS.GoalsExactAchieve indicates the number of objectives for which it is
%      required for the objectives (F) to equal the goals (GOAL). 
%      Such objectives should be partitioned into the first few 
%      elements of F.
%      The remaining parameters determine tolerance settings.
%          
%
%
defaultopt = struct( ...
    'DerivativeCheck','off', ...
    'Diagnostics','off', ...
    'DiffMaxChange',1e-1, ...
    'DiffMinChange',1e-8, ...
    'Display','final', ...
    'FunValCheck','off', ...
    'GoalsExactAchieve',0, ...
    'GradConstr','off', ...
    'GradObj','off', ...
    'Hessian','off', ...
    'LargeScale','off', ...
    'MaxFunEvals','100*numberOfVariables', ...
    'MaxIter',400, ...
    'MaxSQPIter','10*max(numberOfVariables,numberOfInequalities+numberOfBounds)', ...
    'MeritFunction','multiobj', ...
    'NoStopIfFlatInfeas','off', ...
    'OutputFcn',[], ...
    'PhaseOneTotalScaling','off', ...
    'PlotFcns',[], ...
    'RelLineSrchBnd',[], ...
    'RelLineSrchBndDuration',1, ...
    'TolCon',1e-6, ...
    'TolConSQP',1e-6, ...
    'TolFun',1e-6, ...
    'TolX',1e-6, ...
    'TypicalX','ones(numberOfVariables,1)', ...
    'UseParallel','never' ...
    );

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FUN,'defaults')
   x = defaultopt;
   return
end

caller='fgoalattain';

if nargin < 12, options = [];
   if nargin < 11, NONLCON = [];
      if nargin < 10, UB = [];
         if nargin < 9, LB = [];
            if nargin < 8, Beq = [];
               if nargin < 7, Aeq = [];
                  if nargin < 6, B = [];
                     if nargin < 5, A = [];
                        end,end,end,end,end,end,end,end

% Detect problem structure input
problemInput = false;
if nargin == 1
    if isa(FUN,'struct')
        problemInput = true;
        [FUN,x,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error('optim:fgoalattain:InputArg','The input to FGOALATTAIN should be either a structure with valid fields or consist of at least four arguments.');
    end
end

if nargin < 4 && ~problemInput
    error('optim:fgoalattain:NotEnoughInputs','FGOALATTAIN requires four input arguments.')
end

% Check for non-double inputs
% SUPERIORFLOAT errors when superior input is neither single nor double;
% We use try-catch to override SUPERIORFLOAT's error message when input
% data type is integer.
try
    dataType = superiorfloat(x,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB);
    if ~isequal('double', dataType)
        error('optim:fgoalattain:NonDoubleInput', ...
            'FGOALATTAIN only accepts inputs of data type double.')
    end
catch
    error('optim:fgoalattain:NonDoubleInput', ...
        'FGOALATTAIN only accepts inputs of data type double.')
end

initVals.xOrigShape = x;
xnew=[x(:);0];

numberOfVariablesplus1 = length(xnew);
numberOfVariables = numberOfVariablesplus1 - 1;
WEIGHT = WEIGHT(:);
GOAL = GOAL(:);

diagnostics = isequal(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
switch optimget(options,'Display',defaultopt,'fast')
case {'off','none'}
   verbosity = 0;
case 'notify'
   verbosity = 1;  
case 'final'
   verbosity = 2; 
case 'iter'
   verbosity = 3;
otherwise
   verbosity = 2;
end

% Set to column vectors
B = B(:);
Beq = Beq(:);

[xnew(1:numberOfVariables),l,u,msg] = checkbounds(xnew(1:numberOfVariables),LB,UB,numberOfVariables);
if ~isempty(msg)
   EXITFLAG = -2;
   [FVAL,ATTAINFACTOR,LAMBDA] = deal([]);
   OUTPUT.iterations = 0;
   OUTPUT.funcCount = 0;
   OUTPUT.stepsize = [];
   OUTPUT.lssteplength = [];
   OUTPUT.algorithm = 'goal attainment SQP, Quasi-Newton, line_search';
   OUTPUT.firstorderopt = [];
   OUTPUT.constrviolation = [];
   OUTPUT.message = msg;
   x(:) = xnew(1:numberOfVariables);
   if verbosity > 0
      disp(msg)
   end
   return
end

neqgoals = optimget(options, 'GoalsExactAchieve',defaultopt,'fast');
% flags.meritFunction is 1 unless changed by user to fmincon merit function;
% formerly options(7)
% 0 uses the fmincon single-objective merit and Hess; 1 is the default
flags.meritFunction = strcmp(optimget(options,'MeritFunction',defaultopt,'fast'),'multiobj');

lenVarIn = length(varargin);
% goalcon and goalfun also take:
% neqgoals,funfcn,gradfcn,WEIGHT,GOAL,x,errCheck
goalargs = 7; 

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
flags.grad = strcmp(optimget(options,'GradObj',defaultopt,'fast'),'on');
flags.gradconst = strcmp(optimget(options,'GradConstr',defaultopt,'fast'),'on');
flags.hess = strcmp(optimget(options,'Hessian',defaultopt,'fast'),'on');
if flags.hess
   warning('optim:fgoalattain:UserHessNotUsed','FGOALATTAIN does not use user-supplied Hessian.')
   flags.hess = 0;
end

if isempty(NONLCON)
   constflag = 0;
else
   constflag = 1;
end

line_search = strcmp(optimget(options,'LargeScale',defaultopt,'fast'),'off'); % 0 means trust-region, 1 means line-search
if ~line_search
   warning('optim:fgoalattain:NoLargeScale', ...
           'Large-scale algorithm not currently available for this problem type.')
   line_search = 1;
end

% If nonlinear constraints exist, need 
% either both function and constraint gradients, or not
if constflag
   if flags.grad && ~flags.gradconst
      flags.grad = 0;
   elseif ~flags.grad && flags.gradconst
      flags.gradconst = 0;
   end
end

% Convert to inline function as needed
% FUN is called from goalcon; goalfun is based only on x
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
   [funfcn, msg] = optimfcnchk(FUN,'goalcon',length(varargin),funValCheck, ...
       flags.grad,flags.hess);
else
   error('optim:fgoalattain:InvalidFUN', ...
         'FUN must be a function handle or a cell array of two function handles.')
end
% Pass in false for funValCheck argument as goalfun is not a user function
[ffun, msg] = optimfcnchk(@goalfun,'fgoalattain',lenVarIn+goalargs,false,flags.grad);

if constflag % NONLCON is non-empty
   [confcn, msg] = ...
      optimfcnchk(NONLCON,'goalcon',length(varargin),funValCheck,flags.gradconst,0,1);
else
   confcn{1} = '';
end
% Pass in false for funValCheck argument as goalfun is not a user function
[cfun, msg]= optimfcnchk(@goalcon,'fgoalattain',lenVarIn+goalargs,false,flags.gradconst,0,1); 

lenvlb=length(l);
lenvub=length(u);

i=1:lenvlb;
lindex = xnew(i)<l(i);
if any(lindex),
   xnew(lindex)=l(lindex)+1e-4; 
end
i = 1:lenvub;
uindex = xnew(i)>u(i);
if any(uindex)
   xnew(uindex)=u(uindex);
end
x(:) = xnew(1:end-1);
len_user_f = length(GOAL); % Assume the length of GOAL is same as length
                           % of user function; we will verify this later.
if length(WEIGHT) ~= length(GOAL)
     error('optim:fgoalattain:InvalidWeightAndGoalSizes', ...
           'Size of WEIGHT must be equal to the size of GOAL.')
end
 
initVals.g=zeros(numberOfVariablesplus1,1);
initVals.H = [];
errCheck = true; % Perform error checking on initial function evaluations

extravarargin = {neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:}}; 
% Evaluate goal function
switch ffun{1}
    case 'fun'
        initVals.f = ffun{3}(xnew,extravarargin{:});
    case 'fungrad'
        [initVals.f,initVals.g] = ffun{3}(xnew,extravarargin{:});
    otherwise
        error('optim:fgoalattain:InvalidCalltype','Undefined calltype in FGOALATTAIN.')
end


% Evaluate goal constraints
switch cfun{1}
    case 'fun'
        [ctmp,ceqtmp] = cfun{3}(xnew,extravarargin{:});
        initVals.ncineq = ctmp(:);
        initVals.nceq = ceqtmp(:);
        initVals.gnc = zeros(numberOfVariablesplus1,length(initVals.ncineq));
        initVals.gnceq = zeros(numberOfVariablesplus1,length(initVals.nceq));
    case 'fungrad'
        [ctmp,ceqtmp,initVals.gnc,initVals.gnceq] = cfun{3}(xnew,extravarargin{:});
        initVals.ncineq = ctmp(:);
        initVals.nceq = ceqtmp(:);
    otherwise
        error('optim:fgoalattain:InvalidCalltype','Undefined calltype in FGOALATTAIN.')
end

non_eq = length(initVals.nceq);
non_ineq = length(initVals.ncineq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A);

if ~isempty(Aeq) && Aeqcol ~= numberOfVariables
   error('optim:fgoalattain:InvalidSizeOfAeq','Aeq has the wrong number of columns.')
end
if ~isempty(A) && Acol ~= numberOfVariables
   error('optim:fgoalattain:InvalidSizeOfA','A has the wrong number of columns.')
end

just_user_constraints = non_ineq - len_user_f - neqgoals;
OUTPUT.algorithm = 'goal attainment SQP, Quasi-Newton, line_search';  % override nlconst output

if diagnostics > 0
   % Do diagnostics on information so far
   msg = diagnose('fgoalattain',OUTPUT,flags.grad,flags.hess,constflag,flags.gradconst,...
      line_search,options,defaultopt,xnew(1:end-1),non_eq,...
      just_user_constraints,lin_eq,lin_ineq,LB,UB,funfcn,confcn,initVals.f,initVals.g,initVals.H, ...
      initVals.ncineq(1:just_user_constraints),initVals.nceq,initVals.gnc(1:just_user_constraints,:),initVals.gnceq);
end

% Add extra column to account for extra xnew component
A = [A,zeros(lin_ineq,1)];
Aeq = [Aeq,zeros(lin_eq,1)];

% Only need to perform error checking on initial function evaluations
errCheck = false;

% Convert function handles to anonymous functions with additional arguments
% in its workspace.
ffun{3} = @(y,varargin) ffun{3}(y,neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:});
cfun{3} = @(y,varargin) cfun{3}(y,neqgoals,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin{:});

[xnew,ATTAINFACTOR,lambda,EXITFLAG,OUTPUT]=...
   nlconst(ffun,xnew,l,u,full(A),B,full(Aeq),Beq,cfun,options,defaultopt, ...
   verbosity,flags,initVals,varargin{:});

if ~isempty(lambda)
    just_user_constraints = length(lambda.ineqnonlin) - len_user_f - neqgoals;
    lambda.ineqnonlin = lambda.ineqnonlin(1:just_user_constraints);
end
LAMBDA=lambda;

OUTPUT.algorithm = 'goal attainment SQP, Quasi-Newton, line_search';  % override nlconst output

% Evaluate user objective functions
x(:)=xnew(1:end-1);
FVAL = funfcn{3}(x,varargin{:});
