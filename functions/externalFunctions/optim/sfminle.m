function[x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = sfminle(funfcn,x,A,b,verb,options,defaultopt,...
    computeLambda,initialf,initialGRAD,initialHESS,Hstr,varargin)
%SFMINLE Nonlinear minimization with linear equalities.
%
% Locate a local minimizer to 
%
%               min { f(x) :  Ax = b}.
%
%	where f(x) maps n-vectors to scalars.
%

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/12/10 21:50:11 $

%   Initialization
xcurr = x(:); % x has "the" shape; xcurr is a vector
numFunEvals = 1;  % done in calling function fmincon

n = length(xcurr); iter= 0;
header = sprintf(['\n                                Norm of      First-order \n',...
        ' Iteration        f(x)          step          optimality   CG-iterations']);
formatstrFirstIter = ' %5.0f      %13.6g                  %12.3g                ';
formatstr = ' %5.0f      %13.6g  %13.6g   %12.3g     %7.0f';

if n == 0, 
    error('optim:sfminle:InvalidN','n must be positive.')
end
[mm,nn] = size(A);
if n ~= nn
    error('optim:sfminle:AeqAndXInconsistent', ...
          'Column dimension of Aeq is inconsistent with length of X.')
end
m = length(b);
if m ~= mm
    error('optim:sfminle:AeqAndBeqInconsistent', ...
          'Row dimension of Aeq is inconsistent with length of beq.')
end

numberOfVariables = n;

% get options out
pcmtx = optimget(options,'Preconditioner',@preaug) ; %not default yet, so use slow optimget
mtxmpy = optimget(options,'HessMult',defaultopt,'fast');
% Use internal Hessian-multiply function if user does not provide HessMult function 
% or options.Hessian is off
if isempty(mtxmpy) || (~strcmpi(funfcn{1},'fungradhess') && ~strcmpi(funfcn{1},'fun_then_grad_then_hess'))
    mtxmpy = @hmult;
end
pcflags = optimget(options,'PrecondBandWidth',defaultopt,'fast') ;
tol2 = optimget(options,'TolX',defaultopt,'fast') ;
tol1 = optimget(options,'TolFun',defaultopt,'fast') ;
maxiter = optimget(options,'MaxIter',defaultopt,'fast') ;
maxfunevals = optimget(options,'MaxFunEvals',defaultopt,'fast') ;
pcgtol = optimget(options,'TolPCG',defaultopt,'fast') ;  % pcgtol = .1;
pcgtol = min(pcgtol, 1e-2);  % better default
kmax = optimget(options,'MaxPCGIter',defaultopt,'fast') ;
typx = optimget(options,'TypicalX',defaultopt,'fast') ;
if ischar(typx)
    if isequal(lower(typx),'ones(numberofvariables,1)')
        typx = ones(numberOfVariables,1);
    else
        error('optim:sfminle:InvalidTypicalX', ...
              'Option ''TypicalX'' must be a matrix (not a string) if not the default.')
    end
end
DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');

if ischar(kmax)
    if isequal(lower(kmax),'max(1,floor(numberofvariables/2))')
        kmax = max(1,floor(numberOfVariables/2));
    else
        error('optim:sfminle:InvalidMaxPCGIter', ...
              'Option ''MaxPCGIter'' must be an integer value if not the default.')
    end
end
if ischar(maxfunevals)
    if isequal(lower(maxfunevals),'100*numberofvariables')
        maxfunevals = 100*numberOfVariables;
    else
        error('optim:sfminle:InvalidMaxFunEvals', ...
              'Option ''MaxFunEvals'' must be an integer value if not the default.')
    end
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
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

dnewt = []; 
snod = [];
ex = 0; posdef = 1; npcg = 0;
pcgit = 0; delta = 100;nrmsx = 1; ratio = 0; degen = inf; 
dv = ones(n,1);
oval = inf;  gradf = zeros(n,1); newgrad = gradf; Z = []; 

% First-derivative check
if DerivativeCheck
    %
    % Calculate finite-difference gradient
    %
    gf=finitedifferences(xcurr,x,funfcn,[],[],[],initialf,[],...
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

%   Remove (numerical) linear dependencies from A
AA = A; bb = b;
[A,b] = dep(AA,[],bb);
m = size(A,1); mm = size(AA,1);
if verb > 2 && m ~= mm
    fprintf('Linear dependencies detected: %i linear constraint(s) removed.\n',mm-m);
end

%   Get feasible: nearest feas. pt. to xstart
xcurr = feasibl(A,b,xcurr);


% Make xcurr conform to the user's input x
x(:) = xcurr;

if ~isempty(Hstr) % use sparse finite differencing
    
    switch funfcn{1}
        case 'fun'
            error('optim:sfminle:ShouldNotReachThis','should not reach this')
        case 'fungrad'
            val = initialf; gradf(:) = initialGRAD;
        case 'fun_then_grad'
            val = initialf; gradf(:) = initialGRAD;
        otherwise
            error('optim:sfminle:UndefinedCalltypeInFMINCON', ...
                  'Undefined calltype in FMINCON.')
    end
    
    %      Determine coloring/grouping for sparse finite-differencing
    p = colamd(Hstr)'; p = (n+1)*ones(n,1)-p; group = color(Hstr,p);
    H = sfd(x,gradf,Hstr,group,[],DiffMinChange,DiffMaxChange, ...
            funfcn,varargin{:});
    
else % user-supplied computation of H or dnewt
    switch funfcn{1}
        case 'fungradhess'
            val = initialf; gradf(:) = initialGRAD; H = initialHESS;
        case 'fun_then_grad_then_hess'
            val = initialf; gradf(:) = initialGRAD; H = initialHESS;
        otherwise
            error('optim:sfminle:UndefinedCalltypeInFMINCON', ...
                  'Undefined calltype in FMINCON.')
    end
end
[nn,pnewt] = size(gradf);

%   Extract the Newton direction?
if pnewt == 2, dnewt = gradf(1:n,2); end
PT = findp(A);
[g,LZ,UZ,pcolZ,PZ] = project(A,-gradf(1:n,1),PT);
gnrm = norm(g,inf);

if verb > 2
    disp(header)
end

% Initialize the output and plot functions.
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

%   MAIN LOOP: GENERATE FEAS. SEQ.  xcurr(iter) S.T. f(xcurr(iter)) IS DECREASING.
while ~ex
    if ~isfinite(val) || any(~isfinite(gradf))
        error('optim:sfminle:InvalidUserFunction', ...
              '%s cannot continue: user function is returning Inf or NaN values.',funfcn{2})
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
    
    % OutputFcn and PlotFcns call
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
    if gnrm < tol1 && posdef == 1
        ex = 3; EXITFLAG = 1;
        outMessage = ...
          sprintf(['Optimization terminated: first-order optimality less than OPTIONS.TolFun,\n' ...
                   ' and no negative/zero curvature detected in trust region model.']);
        if verb > 1
            disp(outMessage)
        end
    elseif (nrmsx < .9*delta) && (ratio > .25) && (diff < tol1*(1+abs(oval)))
        ex = 1; EXITFLAG = 3;
        outMessage = sprintf(['Optimization terminated:' ...
             ' relative function value changing by less than OPTIONS.TolFun.']);        
        if verb > 1
            disp(outMessage)
        end
        
    elseif iter > 1 && nrmsx < tol2
        ex = 2; EXITFLAG = 2;
        outMessage = sprintf(['Optimization terminated:' ...
                    ' norm of the current step is less than OPTIONS.TolX.']);
        if verb > 1
            disp(outMessage)
        end
        
    elseif iter > maxfunevals
        ex=4; EXITFLAG = 0;
        outMessage = sprintf(['Maximum number of function evaluations exceeded;\n' ...
                            '   increase options.MaxFunEvals.']);
        if verb > 0
            disp(outMessage)
        end
    elseif iter > maxiter
        ex=4; EXITFLAG = 0;
        outMessage = sprintf(['Maximum number of iterations exceeded;\n' ...
                            '   increase options.MaxIter.']);
        if verb > 0
            disp(outMessage)
        end
    end
    
 
    %     Step computation
    if ~ex
        % OutputFcn call
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
        sx = zeros(n,1); 
        oposdef = posdef;
        [sx,snod,qp,posdef,pcgit,Z] = trdg(xcurr,gradf(:,1),H,...
            delta,g,mtxmpy,pcmtx,pcflags,...
            pcgtol,kmax,A,zeros(m,1),Z,dnewt,options,defaultopt,...
            PT,LZ,UZ,pcolZ,PZ,varargin{:});
        
        if isempty(posdef), 
            posdef = oposdef; 
        end
        nrmsx=norm(snod); 
        npcg=npcg + pcgit;
        newx=xcurr + sx; 
        
        % Make newx conform to user's input x
        x(:) = newx;
        %       Evaluate f,g,  and H
        if ~isempty(Hstr) % use sparse finite differencing
            switch funfcn{1}
                case 'fun'
                    error('optim:sfminle:ShouldNotReachThis','should not reach this')
                case 'fungrad'
                    [newval,newgrad(:)] = feval(funfcn{3},x,varargin{:});
                case 'fun_then_grad'
                    newval = feval(funfcn{3},x,varargin{:}); 
                    newgrad(:) = feval(funfcn{4},x,varargin{:});
                otherwise
                    error('optim:sfminle:UndefinedCalltypeInFMINCON', ...
                          'Undefined calltype in FMINCON.')
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
                    error('optim:sfminle:UndefinedCalltypeInFMINCON', ...
                          'Undefined calltype in FMINCON.')
            end
            
        end
        numFunEvals = numFunEvals + 1;
        
        [nn,pnewt] = size(newgrad);
        if pnewt == 2, 
            dnewt = newgrad(1:n,2); 
        end
        aug = .5*snod'*((dv.*abs(newgrad(:,1))).*snod);
        ratio = (newval + aug -val)/qp;
        if (ratio >= .75) && (nrmsx >= .9*delta),
            delta = 2*delta;
        elseif ratio <= .25, 
            delta = min(nrmsx/4,delta/4);
        end
        if newval == inf; 
            delta = min(nrmsx/20,delta/20);
        end
        
        %       Update
        if newval < val
            xcurr = newx; 
            val = newval; 
            gradf = newgrad; 
            H = newH;
            Z = [];
            if pnewt == 2, 
                dnewt = newgrad(1:n,2); 
            end
            g = project(A,-gradf(:,1),PT,LZ,UZ,pcolZ,PZ);
            gnrm = norm(g,inf);
            
            %          Extract the Newton direction?
            if pnewt == 2, 
                dnewt = newgrad(1:n,2); 
            end
        end
        iter = iter+1;
        
    end % if ~ex
end % while

if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,xcurr,xOutputfcn,'done',iter,numFunEvals, ...
        val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,varargin{:});
end


HESSIAN = H;
GRAD = g;
FVAL = val;
if computeLambda
    % Disable the warnings about conditioning for singular and
    % nearly singular matrices
    warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
    warningstate2 = warning('off', 'MATLAB:singularMatrix');

    LAMBDA.eqlin = -A'\gradf;

    % Restore the warning states to their original settings
    warning(warningstate1)
    warning(warningstate2)
    LAMBDA.ineqlin = []; LAMBDA.upper = []; LAMBDA.lower = [];
else
    LAMBDA = [];
end
OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.cgiterations = npcg;
OUTPUT.firstorderopt = gnrm;
OUTPUT.algorithm = 'large-scale: projected trust-region Newton';
OUTPUT.message = outMessage;
x(:) = xcurr;
%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,state,iter,numFunEvals, ...
    val,nrmsx,g,gnrm,pcgit,posdef,ratio,degen,delta,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then
% calls the outputfcn/plotfcns.  
%
% state - can have the values 'init','iter','interrupt', or 'done'. 

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
            error('optim:sfminle:UnknownStateInCALLOUTPUTANDPLOTFCNS', ...
                'Unknown state in CALLOUTPUTANDPLOTFCNS.')
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
            error('optim:sfminle:UnknownStateInCALLOUTPUTANDPLOTFCNS', ...
                'Unknown state in CALLOUTPUTANDPLOTFCNS.')
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

