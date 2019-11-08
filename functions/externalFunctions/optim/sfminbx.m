function [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = sfminbx(funfcn,x,l,u,verb,options,defaultopt,...
    computeLambda,initialf,initialGRAD,initialHESS,Hstr,varargin)
%SFMINBX Nonlinear minimization with box constraints.
%
% Locate a local minimizer to 
%
%               min { f(x) :  l <= x <= u}.
%
%	where f(x) maps n-vectors to scalars.

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/12/10 21:50:10 $

%   Initialization
xcurr = x(:);  % x has "the" shape; xcurr is a vector
n = length(xcurr); 
iter = 0; 
numFunEvals = 1;  % feval of user function in fmincon or fminunc

header = sprintf(['\n                                Norm of      First-order \n',...
        ' Iteration        f(x)          step          optimality   CG-iterations']);
formatstrFirstIter = ' %5.0f      %13.6g                  %12.3g                ';
formatstr = ' %5.0f      %13.6g  %13.6g   %12.3g     %7.0f';
if n == 0, 
    error('optim:sfminbx:InvalidN','n must be positive.')
end
% If called from fminunc, need to fill in l and u.
if isempty(l)
    l = -inf(n,1);
end
if isempty(u)
    u = inf(n,1);
end
arg = (u >= 1e10); arg2 = (l <= -1e10);
u(arg) = inf;
l(arg2) = -inf;
if any(u == l) 
    error('optim:sfminbx:InvalidBounds', ...
          ['Equal upper and lower bounds not permitted in this large-scale method.\n' ...
           'Use equality constraints and the medium-scale method instead.'])
elseif min(u-l) <= 0
    error('optim:sfminbx:InconsistentBounds','Inconsistent bounds.')
end

% get options out
typx = optimget(options,'TypicalX',defaultopt,'fast') ;
% In case the defaults were gathered from calling: optimset('quadprog'):
numberOfVariables = n;
if ischar(typx)
    if isequal(lower(typx),'ones(numberofvariables,1)')
        typx = ones(numberOfVariables,1);
    else
        error('optim:sfminbx:InvalidTypicalX', ...
              'Option ''TypicalX'' must be a matrix (not a string) if not the default.')
    end
end
DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');

% Will be user-settable later:
pcmtx = optimget(options,'Preconditioner',@hprecon) ; % not a default yet

mtxmpy = optimget(options,'HessMult',defaultopt,'fast');
% Use internal Hessian-multiply function if user does not provide HessMult function 
% or options.Hessian is off
if isempty(mtxmpy) || (~strcmpi(funfcn{1},'fungradhess') && ~strcmpi(funfcn{1},'fun_then_grad_then_hess'))
    mtxmpy = @hmult; % to detect name clash with user hmult.m, need this
end

% Handle the output
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot functions
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcn which is needed to support cell array syntax for PlotFcn.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

active_tol = optimget(options,'ActiveConstrTol',sqrt(eps)); % not a default yet, so use slow optimget
pcflags = optimget(options,'PrecondBandWidth',defaultopt,'fast') ;
tol2 = optimget(options,'TolX',defaultopt,'fast') ;
tol1 = optimget(options,'TolFun',defaultopt,'fast') ;
maxiter = optimget(options,'MaxIter',defaultopt,'fast') ;
maxfunevals = optimget(options,'MaxFunEvals',defaultopt,'fast') ;
pcgtol = optimget(options,'TolPCG',defaultopt,'fast') ;  % pcgtol = .1;
kmax = optimget(options,'MaxPCGIter', defaultopt,'fast') ;
if ischar(kmax)
    if isequal(lower(kmax),'max(1,floor(numberofvariables/2))')
        kmax = max(1,floor(numberOfVariables/2));
    else
        error('optim:sfminbx:InvalidMaxPCGIter', ...
              'Option ''MaxPCGIter'' must be an integer value if not the default.')
    end
end
if ischar(maxfunevals)
    if isequal(lower(maxfunevals),'100*numberofvariables')
        maxfunevals = 100*numberOfVariables;
    else
        error('optim:sfminbx:InvalidMaxFunEvals', ...
              'Option ''MaxFunEvals'' must be an integer value if not the default.')
    end
end

dnewt = [];  
ex = 0; posdef = 1; npcg = 0; 

pcgit = 0; delta = 10;nrmsx = 1; ratio = 0;
oval = inf;  g = zeros(n,1); newgrad = g; Z = []; 

% First-derivative check
if DerivativeCheck
    %
    % Calculate finite-difference gradient
    % 
    gf=finitedifferences(xcurr,x,funfcn,[],l,u,initialf,[],...
        [],[],DiffMinChange,DiffMaxChange,typx,[],'all',...
        [],[],[],[],[],[],false,[],varargin{:});
    %
    % Compare user-supplied versus finite-difference gradients
    %
    disp('Function derivative')
    if isa(funfcn{4},'inline')
        graderr(gf, initialGRAD, formula(funfcn{4}));
    else 
        graderr(gf, initialGRAD, funfcn{4});
    end
    numFunEvals = numFunEvals + numberOfVariables;
end

if ~isempty(Hstr)  % use sparse finite differencing
    switch funfcn{1}
        case 'fun'
            error('optim:sfminbx:ShouldNotReachThis','should not reach this')
        case 'fungrad'
            val = initialf; g(:) = initialGRAD;
        case 'fun_then_grad'
            val = initialf; g(:) = initialGRAD;
        otherwise
            if isequal(funfcn{2},'fmincon')
                error('optim:sfminbx:UndefinedFMINCONCalltype','Undefined calltype in FMINCON.')
            else
                error('optim:sfminbx:UndefinedFMINUNCCalltype','Undefined calltype in FMINUNC.')
            end
    end
    % Determine coloring/grouping for sparse finite-differencing
    p = colamd(Hstr)'; 
    p = (n+1) - p;
    group = color(Hstr,p);
    % pass in the user shaped x
    H = sfd(x,g,Hstr,group,[],DiffMinChange,DiffMaxChange,funfcn,varargin{:});
    %
else % user-supplied computation of H or dnewt
    switch funfcn{1}
        case 'fungradhess'
            val = initialf; g(:) = initialGRAD; H = initialHESS;
        case 'fun_then_grad_then_hess'
            val = initialf; g(:) = initialGRAD; H = initialHESS;
        otherwise
            if isequal(funfcn{2},'fmincon')
                error('optim:sfminbx:UndefinedFMINCONCalltype','Undefined calltype in FMINCON.')
            else
                error('optim:sfminbx:UndefinedFMINUNCCalltype','Undefined calltype in FMINUNC.')
            end
    end
end

delbnd = max(100*norm(xcurr),1);
[nn,pp] = size(g);

%   Extract the Newton direction?
if pp == 2, dnewt = g(1:n,2); end
if verb > 2
    disp(header)
end

% Initialize the output function.
if haveoutputfcn || haveplotfcn
   [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xcurr,xOutputfcn,'init',iter,numFunEvals, ...
        val,[],[],[],pcgit,[],[],[],delta,varargin{:});
    if stop
        [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = cleanUpInterrupt(xOutputfcn,optimValues,npcg);
        if verb > 0
            disp(OUTPUT.message)
        end
        return;
    end
end

%   MAIN LOOP: GENERATE FEAS. SEQ.  xcurr(iter) S.T. f(x(iter)) IS DECREASING.
while ~ex
    if ~isfinite(val) || any(~isfinite(g))
        error('optim:sfminbx:InvalidUserFunction', ...
              '%s cannot continue: user function is returning Inf or NaN values.',funfcn{2})
    end
    
    %     Update
    [v,dv] = definev(g(:,1),xcurr,l,u); 
    gopt = v.*g(:,1); gnrm = norm(gopt,inf);
    r = abs(min(u-xcurr,xcurr-l)); degen = min(r + abs(g(:,1)));
    % If no upper and lower bounds (e.g., fminunc), then set degen flag to -1.
    if all( (l == -inf) & (u == inf) )
        degen = -1; 
    end
    
    % Display
    if verb > 2
        if iter > 0
            currOutput = sprintf(formatstr,iter,val,nrmsx,gnrm,pcgit);
        else
            currOutput = sprintf(formatstrFirstIter,iter,val,gnrm);
        end
        disp(currOutput);
    end
    
    % OutputFcn call
    if haveoutputfcn || haveplotfcn
       [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xcurr,xOutputfcn,'iter',iter,numFunEvals, ...
            val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,varargin{:});
        if stop  % Stop per user request.
            [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
                cleanUpInterrupt(xOutputfcn,optimValues,npcg);
            if verb > 0
                disp(OUTPUT.message)
            end
            return;
        end
    end
    
    %     TEST FOR CONVERGENCE
    diff = abs(oval-val); 
    oval = val; 
    if ((gnrm < tol1) && (posdef ==1) )
        ex = 3;
        outMessage = ...
          sprintf(['Optimization terminated: first-order optimality less than OPTIONS.TolFun,\n', ...
                   ' and no negative/zero curvature detected in trust region model.']);        
        if verb > 1
            disp(outMessage)
        end
    elseif (nrmsx < .9*delta) && (ratio > .25) && (diff < tol1*(1+abs(oval)))
        ex = 1;
        outMessage = sprintf(['Optimization terminated:' ... 
                  ' relative function value changing by less than OPTIONS.TolFun.']);
        if verb > 1
            disp(outMessage)
        end
    elseif (iter > 1) && (nrmsx < tol2) 
        ex = 2;
        outMessage = sprintf(['Optimization terminated:' ...
                    ' norm of the current step is less than OPTIONS.TolX.']);
        if verb > 1
            disp(outMessage)
        end
    elseif numFunEvals > maxfunevals
        ex = 4;
        outMessage = sprintf(['Maximum number of function evaluations exceeded;\n' ...
                            '   increase options.MaxFunEvals.']);
        if verb > 0
            disp(outMessage)
        end
    elseif iter > maxiter
        ex = 4;
        outMessage = sprintf(['Maximum number of iterations exceeded;\n' ...
                            '   increase options.MaxIter.']);
        if verb > 0
            disp(outMessage)
        end
    end
      
    %     Step computation
    if ~ex
        if haveoutputfcn % Call output functions (we don't call plot functions with 'interrupt' flag)
            [unused1, unused2, stop] = callOutputAndPlotFcns(outputfcn,{},xcurr,xOutputfcn,'interrupt',iter,numFunEvals, ...
                val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,varargin{:});
            if stop  % Stop per user request.
                [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
                    cleanUpInterrupt(xOutputfcn,optimValues,npcg);
                if verb > 0
                    disp(OUTPUT.message)
                end
                return;
            end
        end
        %       Determine trust region correction
        dd = abs(v); D = sparse(1:n,1:n,full(sqrt(dd))); 
        theta = max(.95,1-gnrm);  
        oposdef = posdef;
        [sx,snod,qp,posdef,pcgit,Z] = trdog(xcurr, g(:,1),H,D,delta,dv,...
            mtxmpy,pcmtx,pcflags,pcgtol,kmax,theta,l,u,Z,dnewt,'hessprecon',varargin{:});
        if isempty(posdef), posdef = oposdef; end
        nrmsx=norm(snod); npcg=npcg + pcgit;
        newx=xcurr + sx;
        
        %       Perturb?
        [pert,newx] = perturb(newx,l,u); 
        
        
        % Make newx conform to user's input x
        x(:) = newx;
        % Evaluate f, g, and H
        if ~isempty(Hstr) % use sparse finite differencing
            switch funfcn{1}
                case 'fun'
                    error('optim:sfminbx:ShouldNotReachThis','should not reach this')
                case 'fungrad'
                    [newval,newgrad(:)] = feval(funfcn{3},x,varargin{:});
                case 'fun_then_grad'
                    newval = feval(funfcn{3},x,varargin{:}); 
                    newgrad(:) = feval(funfcn{4},x,varargin{:});
                otherwise
                    if isequal(funfcn{2},'fmincon')
                        error('optim:sfminbx:UndefinedCalltypeInFMINCON', ...
                              'Undefined calltype in FMINCON.')
                    else
                        error('optim:sfminbx:UndefinedCalltypeInFMINUNC', ...
                              'Undefined calltype in FMINUNC.')
                    end
            end
            newH = sfd(x,newgrad,Hstr,group,[],DiffMinChange,DiffMaxChange, ...
                       funfcn,varargin{:});
            
        else % user-supplied computation of H or dnewt
            switch funfcn{1}
                case 'fungradhess'
                    [newval,newgrad(:),newH] = feval(funfcn{3},x,varargin{:});
                case 'fun_then_grad_then_hess'
                    newval = feval(funfcn{3},x,varargin{:}); 
                    newgrad(:) = feval(funfcn{4},x,varargin{:});
                    newH = feval(funfcn{5},x,varargin{:});
                otherwise
                    if isequal(funfcn{2},'fmincon')
                        error('optim:sfminbx:UndefinedCalltypeInFMINCON', ...
                              'Undefined calltype in FMINCON.')
                    else
                        error('optim:sfminbx:UndefinedCalltypeInFMINUNC', ...
                              'Undefined calltype in FMINUNC.')
                    end
            end
            
        end
        numFunEvals = numFunEvals + 1;
        [nn,pp] = size(newgrad);
        aug = .5*snod'*((dv.*abs(newgrad(:,1))).*snod);
        ratio = (newval + aug -val)/qp;
        
        if (ratio >= .75) && (nrmsx >= .9*delta)
            delta = min(delbnd,2*delta);
        elseif ratio <= .25
            delta = min(nrmsx/4,delta/4);
        end
        if newval == inf
            delta = min(nrmsx/20,delta/20);
        end
        
        %       Update
        if newval < val
            xcurr=newx; val = newval; g= newgrad; H = newH;
            Z = [];
            
            %          Extract the Newton direction?
            if pp == 2, dnewt = newgrad(1:n,2); end
        end
        iter = iter + 1;
    end % if ~ex
end % while

if haveoutputfcn || haveplotfcn
   callOutputAndPlotFcns(outputfcn,plotfcns,xcurr,xOutputfcn,'done',iter,numFunEvals, ...
        val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,varargin{:});
    % Do not check value of 'stop' as we are done with the optimization
    % already.
end
HESSIAN = H;
GRAD = g;
FVAL = val;
LAMBDA = [];

if ex == 3
  EXITFLAG = 1;  % first order optimality conditions hold
elseif ex == 1
  EXITFLAG = 3;  % relative change in f small
elseif ex == 2
  EXITFLAG = 2;  % norm of step small
elseif ex == 4   % MaxIter or MaxFunEval
  EXITFLAG = 0;
else             % ex = 10
  EXITFLAG = -1; % Exiting per request
end

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.cgiterations = npcg;
OUTPUT.firstorderopt = gnrm;
OUTPUT.algorithm = 'large-scale: trust-region reflective Newton'; 
OUTPUT.message = outMessage;

x(:) = xcurr;
if computeLambda
    g = full(g);
    
    LAMBDA.lower = zeros(length(l),1);
    LAMBDA.upper = zeros(length(u),1);
    argl = logical(abs(xcurr-l) < active_tol);
    argu = logical(abs(xcurr-u) < active_tol);
    
    LAMBDA.lower(argl) = (g(argl));
    LAMBDA.upper(argu) = -(g(argu));
    LAMBDA.ineqlin = []; LAMBDA.eqlin = []; LAMBDA.ineqnonlin=[]; LAMBDA.eqnonlin=[];
else
    LAMBDA = [];   
end

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,state, ... 
    iter,numFunEvals, ...
    val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% The input STATE can have the values 'init','iter','interrupt', or 'done'. 
%
% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.

optimValues.iteration = iter;
optimValues.funccount = numFunEvals;
optimValues.fval = val;
optimValues.stepsize = nrmsx;
optimValues.gradient = g;
optimValues.firstorderopt = gnrm;
optimValues.cgiterations = pcgit; 
optimValues.positivedefinite = posdef;
optimValues.ratio = ratio;
optimValues.degenerate = min(degen,1);
optimValues.trustregionradius = delta;
optimValues.procedure = '';

xOutputfcn(:) = xvec;  % Set x to have user expected size
stop = false;
if ~isempty(outputfcn)
    switch state
        case {'iter','init','interrupt'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('optim:sfminbx:UnknownStateInCALLOUTPUTANDPLOTFCNS','Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('optim:sfminbx:UnknownStateInCALLOUTPUTANDPLOTFCNS','Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end
%--------------------------------------------------------------------------
function [x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = cleanUpInterrupt(xOutputfcn,optimValues,npcg)
% CLEANUPINTERRUPT updates or sets all the output arguments of SFMINBX when the optimization 
% is interrupted.  The HESSIAN and LAMBDA are set to [] as they may be in a
% state that is inconsistent with the other values since we are
% interrupting mid-iteration.

x = xOutputfcn; 
FVAL = optimValues.fval;
EXITFLAG = -1; 
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.stepsize = optimValues.stepsize;
OUTPUT.algorithm = 'large-scale: trust-region reflective Newton'; 
OUTPUT.firstorderopt = optimValues.firstorderopt; 
OUTPUT.cgiterations = npcg; % total number of CG iterations
OUTPUT.message = 'Optimization terminated prematurely by user.';
GRAD = optimValues.gradient;
HESSIAN = []; % May be in an inconsistent state
LAMBDA = []; % May be in an inconsistent state


