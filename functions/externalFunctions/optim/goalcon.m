function [GS,GSeq,GGS,GGSeq] = goalcon(V,neqgoal,funfcn,confcn,WEIGHT,GOAL,x,errCheck,varargin)
%GOALCON Utility function to translate gradient in goal-attainment problem.
%   Intermediate function used to translate goal attainment
%   problem into constrained optimization problem.
%   Used by FGOALATTAIN and FMINIMAX.
%
%   See also GOALFUN.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/06/14 05:19:03 $

nx=length(V)-1;
x(:)=V(1:nx);
lambda=V(nx+1);

% Compute the constraints
f = []; gf = [];  % Tell parser that f and g are variables.
gc = []; gceq = []; c = []; ceq = [];

% Determine whether or not to check for errors in objective functions
if ~errCheck         
    switch funfcn{1} % evaluate objective functions and possibly gradients
        case 'fun'
            f = funfcn{3}(x,varargin{:});  
        case 'fungrad'
            [f,gf] = funfcn{3}(x,varargin{:});
        case 'fun_then_grad'
            f = funfcn{3}(x,varargin{:});  
            gf = funfcn{4}(x,varargin{:});
        otherwise
            error('optim:goalcon:UndefinedCalltypeFgoalattain','Undefined calltype in FGOALATTAIN.');
    end

    % Evaluate constraints
    switch confcn{1}
        case 'fun'
            [ctmp,ceqtmp] = confcn{3}(x,varargin{:});
            c = ctmp(:); ceq = ceqtmp(:);
            gc = [];
            gceq = [];
        case 'fungrad'
            [ctmp,ceqtmp,gc,gceq] = confcn{3}(x,varargin{:});
            c = ctmp(:); ceq = ceqtmp(:);
        case 'fun_then_grad'
            [ctmp,ceqtmp] = confcn{3}(x,varargin{:});
            c = ctmp(:); ceq = ceqtmp(:);
            [gc,gceq] = confcn{4}(x,varargin{:});
        case ''
            c = []; ceq = [];
            gc = [];
            gceq = [];
        otherwise
            error('optim:goalcon:UndefinedCalltypeFgoalattain','Undefined calltype in FGOALATTAIN.');
    end

    nGoalCnstr = length(f); % number of goal attainment constraints 
else
    % Evaluate objective functions and check for internal errors
    switch funfcn{1} % evaluate function and possibly gradients
        case 'fun'
            try
                f = funfcn{3}(x,varargin{:});  
            catch
                error('optim:goalcon:ObjectiveError', ...
                ['The optimization cannot continue because user supplied objective function' ...
                ' failed with the following error:\n%s'], lasterr)
            end
        case 'fungrad'
            try
                [f,gf] = funfcn{3}(x,varargin{:});
            catch
                error('optim:goalcon:ObjectiveError', ...
                ['The optimization cannot continue because user supplied objective function' ...
                ' failed with the following error:\n%s'], lasterr)
            end
        case 'fun_then_grad'
            try
                f = funfcn{3}(x,varargin{:});              
            catch
                error('optim:goalcon:ObjectiveError', ...
                ['The optimization cannot continue because user supplied objective function' ...
                ' failed with the following error:\n%s'], lasterr)
            end
            try
                gf = funfcn{4}(x,varargin{:});
            catch
                error('optim:goalcon:ObjectiveError', ...
                ['The optimization cannot continue because user supplied objective gradient' ...
                ' function failed with the following error:\n%s'], lasterr)
            end
        otherwise
            error('optim:goalcon:UndefinedCalltypeFgoalattain','Undefined calltype in FGOALATTAIN.');
    end

    % Evaluate constraints and check for internal errors
    switch confcn{1}
        case 'fun'
            try
                [ctmp,ceqtmp] = confcn{3}(x,varargin{:});
            catch
                error('optim:goalcon:NonlconErrorError', ...
                ['The optimization cannot continue because user supplied nonlinear constraint function' ...
                ' failed with the following error:\n%s'], lasterr)
            end
            c = ctmp(:); ceq = ceqtmp(:);
            gc = [];
            gceq = [];
        case 'fungrad'
            try
                [ctmp,ceqtmp,gc,gceq] = confcn{3}(x,varargin{:});
            catch
                error('optim:goalcon:NonlconErrorError', ...
                ['The optimization cannot continue because user supplied nonlinear constraint function' ...
                ' failed with the following error:\n%s'], lasterr)
            end
            c = ctmp(:); ceq = ceqtmp(:);
        case 'fun_then_grad'
            try
                [ctmp,ceqtmp] = confcn{3}(x,varargin{:});
            catch
                error('optim:goalcon:NonlconErrorError', ...
                ['The optimization cannot continue because user supplied nonlinear constraint function' ...
                ' failed with the following error:\n%s'], lasterr)
            end
            c = ctmp(:); ceq = ceqtmp(:);
            try
                [gc,gceq] = confcn{4}(x,varargin{:});
            catch
                error('optim:goalcon:NonlconErrorError', ...
                ['The optimization cannot continue because user supplied nonlinear constraint gradient' ...
                ' function failed with the following error:\n%s'], lasterr)
            end
        case ''
            c = []; ceq = [];
            gc = [];
            gceq = [];
        otherwise
            error('optim:goalcon:UndefinedCalltypeFgoalattain','Undefined calltype in FGOALATTAIN.');
    end

    nGoalCnstr = length(f); % number of goal attainment constraints 

    % Initial check of user-functions. The errCheck flag is set accordingly
    % in calling functions. Check the size of user-supplied objective
    % function gradients
    if ~isempty(gf) && ~all(size(gf) == [nx,nGoalCnstr])
        error('optim:goalcon:InvalidSizeOfGF', ...
            'Gradient of the objectives expected to be %d-by-%d.',nx,nGoalCnstr)
    end

    % Check that the length of user objective function is equal to length(GOAL)
    if length(GOAL) ~= nGoalCnstr
        error('optim:goalcon:InvalidGoalAndFunSizes', ...
            'Size of GOAL must be equal to the size of F returned by FUN.')
    end
    
    non_eq = length(ceq);
    non_ineq = length(c);
    [cgrow, cgcol]= size(gc);
    [ceqgrow, ceqgcol]= size(gceq);
    
    % Check the size of user-supplied non-linear constraint gradients
    if ~isempty(gc) && (cgrow~=nx || cgcol~=non_ineq)
        error('optim:goalcon:InvalidSizeOfGC', ...
            'Gradient of the nonlinear inequality constraints expected to be %d-by-%d.',nx,non_ineq)
    end
    if ~isempty(gceq) && (ceqgrow~=nx || ceqgcol~=non_eq)
        error('optim:goalcon:InvalidSizeOfGCeq', ... 
            'Gradient of the nonlinear equality constraints expected to be %d-by-%d.',nx,non_eq)
    end
end

% Calculate goal attainment constraints ( f(i)-GOAL(i) )/WEIGHT(i)-lambda <= 0
GS = zeros(nGoalCnstr+neqgoal,1);
for i=1:nGoalCnstr
     if WEIGHT(i)~=0
       diff=f(i)-GOAL(i);
       GS(i)=sign(real(diff))*norm(diff)/WEIGHT(i)-lambda; 
       if i<=neqgoal % neqgoal comes from options.GoalsExactAchieve 
          GS(i+nGoalCnstr)=-GS(i)-2*lambda; % f(i)+lambda*WEIGHT>=GOAL
       end
     else % hard constraint
       GS(i)=f(i)-GOAL(i);
       if i<=neqgoal 
          GS(i+nGoalCnstr)=-GS(i)-1e-10; % f(i)>=GOAL(i)-1e-10
       end
     end
end

% Append goal attainment constraint at the end of inequalities vector
GS=[c(:);GS]; 

% Equalities and gradient matrix of equalities
GSeq = ceq(:);
size_ceq = size(GSeq);
GGSeq = [gceq; zeros(1,size_ceq(1))]; % add zero row (derivatives w.r.t. lambda)

if isempty(gf) && isempty(gc)
   GGS=[]; GGSeq = [];
elseif (isempty(gf) && ~isempty(gc)) 
   error('optim:goalcon:GradFunRequiredForGradCon','Must provide gradient of function if gradient of constraints to be used.');
else % grad of f available, grad of inequalities available or not
   
   % Gradient matrix of inequality constraints: 
   % [grad_x of inequalities c(x)          | grad_x of attainment constraint]
   % [zero row (derivatives w.r.t. lambda) | -1 -1 . . . . . .  . . . . ..-1]
   GL = -ones(1,nGoalCnstr+neqgoal);
   for i=1:nGoalCnstr
      if WEIGHT(i)~=0
         gf(:,i)=gf(:,i)/WEIGHT(i);
         if i<=neqgoal,
            gf(:,i+nGoalCnstr)=-gf(:,i);
         end 
      else % hard constraint
         GL(1,i)=0; 
      end
   end
   
   GGS=[gf;GL];
   
   sizegc=size(gc);
   % Put gc first
   if sizegc(1)>0, 
      GGS=[[gc;zeros(1,sizegc(2))],GGS]; 
   end
end
