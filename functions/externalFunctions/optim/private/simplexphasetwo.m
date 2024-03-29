function [xsol, fval, dualvars, exitflag, niters, basicVarIdx, nonbasicVarIdx] ...
    = simplexphasetwo(c,A,b,lb,ub,basicVarIdx,nonbasicVarIdx,x0,maxiter,...
    tol,verbosity,computeLambda)
%SIMPLEXPHASETWO Phase two of simplex method for linear programming.
%   X = SIMPLEXPHASETWO(c,A,b,lb,ub) solves linear programming problems in
%   general form with upper bound and lower bound:
%              min  c'* x 
%              s.t. A * x = b,  
%                 lb <= x <= ub.
%
%   [X,FVAL] = SIMPLEXPHASETWO(c,A,b,lb,ub) returns the value of the 
%   objective function at X: FVAL = c'*X.
%
%   [X,FVAL,EXITFLAG] = SIMPLEXPHASETWO(c,A,b,lb,ub) returns EXITFLAG that 
%   describes the exit solution condition of SIMPLEXPHASETWO.
%   If EXITFLAG is:
%      1    then SIMPLEXPHASETWO converged with a solution X.
%      0    then the maximum number of iterations was exceeded
%      < 0  then the problem is unbounded, infeasible, or 
%                SIMPLEXPHASETWO failed to converge with a solution X.
%      -1   then the problem is infeasible
%      -2   then the problem is unbounded
%      -3   then the problem is degenerate
%
%   [X,FVAL,EXITFLAG, NITERS] = SIMPLEXPHASETWO(c,A,b) returns a number of 
%   total iterations NITERS used in the simplex method.
%
%   Input:  A  ... constraint matrix
%           b  ... the right hand side of the constraints
%           c  ... the objective function (coefficient vector)
%           lb ... the lower bound on the variable x
%           ub ... the upper bound on the variable x
%
%   Output: xsol     ... optimal solution
%           fval     ... optimal value
%           exitflag ... the status of the solution

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2007/10/10 20:53:40 $

% -------Set up a feasible basis solution (starting vertex)-------------
% Set up the input problem data: c, A, b, lb, ub and
%        basicVarIdx, nonbasicVarIdx; B, N and c_B, c_N; ub_B, lb_B; and ub_N, lb_N.
if nargin < 11
    verbosity = 0;
end

[m, n] = size(A);

if n > 0 
    if isempty(lb) 
        lb = -inf(n,1);
    end
    if isempty(ub)
        ub = inf(n,1);
    end
else
    error('optim:simplexphasetwo:EmptyA', ...
          'Exiting form SIMPLEXPHASETWO: the constraint matrix A is empty.');
end

c  = c(:);
b  = b(:);
lb = lb(:);
ub = ub(:);
x0 = x0(:);
basicVarIdx = basicVarIdx(:);
nonbasicVarIdx = nonbasicVarIdx(:);

if ( nnz(basicVarIdx) > 0 ) 
    Basis    = sparse(A(:, basicVarIdx)); % Converting the Basis to sparse.
    Nonbasis = A(:, nonbasicVarIdx);
    
    x_B = x0(basicVarIdx);
    x_N = x0(nonbasicVarIdx);
    
    c_B = c(basicVarIdx);
    c_N = c(nonbasicVarIdx);
    
    ub_B = ub(basicVarIdx);
    ub_N = ub(nonbasicVarIdx);
    lb_B = lb(basicVarIdx);
    lb_N = lb(nonbasicVarIdx);    
end

tol2 = tol * 1.e-2;
if verbosity >= 5
    % Check the output data
    ubcheck = max ( [x_B - ub_B; x_N - ub_N] );
    lbcheck = max ( [lb_B - x_B; lb_N - x_N] );
    if (ubcheck > tol2) || (lbcheck > tol2)
        error('optim:simplexphasetwo:BndryCondViolated', ...
              'Boundary condition is violated by %e.', max(ubcheck, lbcheck) );
    end
end

% Initialize the output variables, display the starting iteration
dualvars = struct('y', [], 'z', [], 'w', []);
niters = 0;
xsol = [x_B; x_N];
fval = c_B' * x_B + c_N' * x_N;
if verbosity >= 2
    disp( sprintf('\nPhase 2: Minimize using simplex.') );
    disp( sprintf('      Iter            Objective              Dual Infeasibility ') );
    disp( sprintf('                        f''*x                   A''*y+z-w-f') );
end

if verbosity >= 5
    disp( sprintf('    [basicVarIdx      c_B      x_B      lb_B      ub_B] = ') );
    disp(              full([basicVarIdx      c_B      x_B      lb_B      ub_B]) );
    disp( sprintf('    [nonbasicVarIdx      c_N      x_N      lb_N      ub_N] = ') );
    disp(              full([nonbasicVarIdx      c_N      x_N      lb_N      ub_N]) );
    disp( sprintf('*********************************') );
end

% Test optimality condition and update Basis solution if necessary----------

% Initialization of constants for exitflag
Converged   = 1;
ExcdMaxiter = 0;
%Infeasible  = -1;
Unbounded   = -2; % Will be remapped to -3 in linprog.m
Degenerate  = -3;
Unset       = inf;

% Initialization for the while loop
exitflag = Unset;
tol2 = tol^2;

% sameBasis=true signals whether the entering variable switches between its bounds
sameBasis = false;  

% Disable the warnings about conditioning for singular and
% nearly singular matrices
warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
warningstate2 = warning('off', 'MATLAB:singularMatrix');

% Initialization for Updates
RefactorNow = true;
calculateReducedCost = true;
% eiT is a sparse row-vector which will contain a single "1" entry 
eiT = sparse(1,m);

% Initialization for Devex Pricing
frameWorkIdx = nonbasicVarIdx; % length(nonbasicVarIdx) = n-m
sqrweights = ones(n-m,1); % initialize the square of weights.
% allIndexEnter is a vector indicating the positions of the entering variables 
% in the nonbasic variable stack.
allIndexEnter = sparse(n-m,1); 
% basicFrameWorkVars is a vector indicating which variables in the current
% framework are basic.
basicFrameWorkVars = sparse(m,1);

while (niters < maxiter)
    % Solve the system y * Basis = c_B', 
    if ~sameBasis
        if RefactorNow        
            [L,U,P,Q] = lu(Basis);           
            RefactorNow = false;
            % ETA is a matrix that is used to update the Basis from
            % one iteration to the next. It is an identity matrix with its
            % single column replaced by z(indexLeave). The structure eta
            % contains the information necessary to construct this ETA
            % matrices between refactorizations.            
            eta.cind = []; % indices of the replaced columns.
            eta.matrix = []; % the collection of replaced columns.
            eta.pivot = []; % the collection of pivot elements in the replaced columns.
        end
        if calculateReducedCost
            c_BT = sparse(c_B');
            if ~isempty(eta.cind)
                c_BT = backeta(eta.cind,eta.pivot,eta.matrix,c_BT);
            end
            y = (((c_BT*Q)/U)/L)*P;
            y(abs(y)<1e-10) = 0; 
            dualf = c_N - (y * Nonbasis)';
            calculateReducedCost = false;
        else % Update Reduced Cost
            % Details in I. Maros (2003), Computational Techniques of the
            % Simplex Method, p. 183, Section 9.4.3.
            eiT(indexLeave) = 1;
            tmp_eiT = eiT;
            if ~isempty(eta.cind)
                tmp_eiT = backeta(eta.cind,eta.pivot,eta.matrix,eiT);
            end
            rhop = (((tmp_eiT*Q)/U)/L)*P;
            rhop(abs(rhop)<1e-10) = 0; 
            alpha_N = (rhop * Nonbasis)';
            dualf_indexEnter = dualf(indexEnter);
            dualf = dualf - dualf_indexEnter*alpha_N;
            dualf(indexEnter) = - dualf_indexEnter/z(indexLeave);
            % Prepare eiT for the next iteration
            eiT(indexLeave) = 0;

            % Updating the weights for Devex Pricing
            % Details in P. M. J. Harris (1973) "Pivot Selection
            % Methods of the Devex LP Code, Math. Prog. v. 5, pp. 1-28.
            
            % Note that the update of the weights can be performed as soon
            % as indexEnter and indexLeave is known at the end of an
            % iteration. However to keep the flow of the algorithm the same
            % as the original, the update is performed at the beginning of
            % the following iteration.
            
            % Keeping track of the basic variables that are in the framework
            allIndexEnter(indexEnter) = 1;
            % The variable entering the basis will enter at the location of the
            % leaving variable. If this new basic variable is in the
            % current framework make the corresponding location in
            % basicFrameWorkVars nonzero, otherwise the entering variable
            % (at location indexLeave) is not in the current framework and
            % its entry in basicFrameWorkVars should be zero.
            if ismember(basicVarIdx(indexLeave),frameWorkIdx(allIndexEnter>0))
                basicFrameWorkVars(indexLeave) = indexLeave;
            else
                basicFrameWorkVars(indexLeave) = 0;
            end
            % Calculating the weight for the column entering the basis
            % exactly
            sqrweightExact = sum(z(basicFrameWorkVars>0).^2);

            % If sqrweightExact is 1/4th or less of the sqrweights(indexEnter)
            if sqrweights(indexEnter)/sqrweightExact > 4
                % Resetting the framework
                frameWorkIdx = nonbasicVarIdx;
                sqrweights = ones(n-m,1);
                allIndexEnter(allIndexEnter>0) = 0;
                basicFrameWorkVars(basicFrameWorkVars>0) = 0;
            else
                sqrweights = max(sqrweights,sqrweightExact.*alpha_N.^2);
                % The location of the entering variable is occupied by the leaving
                % variable in the nonbasicvarIdx
                sqrweights(indexEnter) = max(1,sqrweightExact/z(indexLeave)^2); % Note that abs(z(indexLeave)) > 1e-8 holds.

            end            
        end
    end
    
    % Choose the entering variable and check the optimality condition    
    
    enteringCandidates = find( ((dualf >= tol) & (x_N - lb_N > tol2)) | ((dualf <= -tol) & (tol2 < ub_N - x_N)) );
    
    if verbosity >= 2
        disp( sprintf('%8d         %12.6g                %12.6g', niters, full(fval), norm(dualf(enteringCandidates))) );
    end
    
    if isempty( enteringCandidates ) % No entering variable exists with updated dualf,
        % double check with calculated (using the latest factorizations L
        % and U) dualf.
        % Note that, y can be calculated from scratch as y = c_B'/Basis.
        c_BT = sparse(c_B');
        if ~isempty(eta.cind)
            c_BT = backeta(eta.cind,eta.pivot,eta.matrix,c_BT);
        end
        y = (((c_BT*Q)/U)/L)*P;
        y(abs(y)<1e-10) = 0; 
        dualf = c_N - (y * Nonbasis)';
        enteringCandidates = find( ((dualf >= tol) & (x_N - lb_N > tol2)) | ...
            ((dualf <= -tol) & (tol2 < ub_N - x_N)) );
    end
    
    if isempty( enteringCandidates ) % No entering variable exists, already optimal                                                       
        exitflag = Converged;
        % Final output information
        if verbosity >= 4
            disp('  Converged to the optimal solution.');  
        end
        % Sort the optimal solution according to the original problem
        [tmp, order] = sort([basicVarIdx; nonbasicVarIdx]);
        xsol = xsol(order);
        
        % Dual variables for lambda: y, zdual, wdual
        if computeLambda == 1
            y = y(:); % Corresponding to both inequalities and equalities
            dualvars.y = y;
            
            tol3 = max(10e-6,tol);
            indgtlo = (xsol > lb + tol3);
            nnzgtlo = nnz(indgtlo);
            zdual = zeros(size(indgtlo)); % Set the size of zdual right, catch the zeros.
            zdual(indgtlo,1) = zeros(nnzgtlo, 1);
            indeqlo = ( abs(x_N -lb_N) < tol3 );
            zdual(nonbasicVarIdx(indeqlo), 1) = dualf(indeqlo);
            dualvars.z = zdual;
            
            indltup = (xsol < ub - tol3); 
            nnzltup = nnz(indltup); 
            wdual = zeros(size(indltup)); % Make assignment to zeros that the size is right
            wdual(indltup,1) = zeros(nnzltup,1);
            indequp = ( abs(x_N - ub_N) < tol3 );
            wdual(nonbasicVarIdx(indequp), 1) = -dualf(indequp);
            dualvars.w = wdual;
            
            if verbosity >= 4
                disp( sprintf('  The norm of the dual feasibility norm(A''y-w +z -c) = %8.2e', ...
                    norm(A'*y - wdual + zdual - c) ) );
            end
        end
        % Restore the warning states to their original settings
        warning(warningstate1)
        warning(warningstate2)
        return; 
    end
    
    % indexEnter is the chosen index of the entering variable 
    % Devex pricing - Details in P. M. J. Harris (1973) "Pivot Selection
    % Methods of the Devex LP Code, Math. Prog. v. 5, pp. 1-28.
    [vl,ind_enteringCandidates] = max((dualf(enteringCandidates).^2)./sqrweights(enteringCandidates));
    indexEnter = enteringCandidates(ind_enteringCandidates);
    
    
    niters = niters + 1;
    
    % To avoid small values (i.e. < 1e-8) on the diagonal of the eta matrix a
    % loop is added here to change indexEnter until there is no entering
    % candidate. What will happen then?
    acceptCandidate = false;
    
    while ~acceptCandidate

        % Solve the system Basis * z = N_{.,k}
        N_indexEnter = sparse(Nonbasis(:, indexEnter)); % 
        z = Q*(U\(L\(P*N_indexEnter)));
        if ~isempty(eta.cind)
            z = forwardeta(eta.cind,eta.pivot,eta.matrix,z);
        end

        % Choose the leaving variable and update the solution
        indexLeave = 0;
        minUb      = inf;
        zind       = find( abs(z) > tol);
        indnegz    = z(zind) < - tol;
        % tolprt: Perturbation amount for Harris' test        
        tolprt       = max(10e-8,min(10e-6,tol));
        
        if ( dualf(indexEnter) <= -tol)  && ( tol2 < ub_N(indexEnter) - x_N(indexEnter) )
            
            ub_delta = ub_N(indexEnter) - x_N(indexEnter);

            % Find index that gives the tightest upper bound
            if ~isempty(zind)
                
                % Harris's test: Pass 1
                bounds = lb_B(zind) - tolprt;
                bounds(indnegz) = ub_B( zind(indnegz) ) + tolprt;
                tB = (x_B(zind) - bounds)./ z(zind);
                minUb = min(tB);
                
                % Harris's test: Pass 2 
                if minUb <= ub_delta
                    bounds = lb_B(zind);
                    bounds(indnegz) = ub_B( zind(indnegz) );
                    tB = (x_B(zind) - bounds)./ z(zind);
                    lvcands = (tB <= minUb);
                    % Find the index of the variables with the maximum pivot
                    % value in column z among the candidates.
                    [pivotMax_unused,lv] = max(abs(z(zind(lvcands))));
                    % Get the first lv indices corresponding to the nonzero
                    % entries of lvcands.
                    indlvcands = find(lvcands,lv);
                    % The final index corresponds to the index of the leaving
                    % variable.
                    lv = indlvcands(end);
                    minUb = max(0,tB(lv));
                    indexLeave = zind(lv);
                end
            end

            if minUb < ub_delta
                delta = minUb;
                sameBasis = false;
            else
                delta = ub_delta;
                sameBasis = true;
            end

            % Update on solution x_N and x_B
            if delta < inf && delta > tol2
                x_N(indexEnter) = x_N(indexEnter) + delta;
                x_B = x_B - delta * z;
            end

        elseif ( dualf(indexEnter) >= tol )  &&  ( x_N(indexEnter) - lb_N(indexEnter) > tol2 )
            ub_delta =  x_N(indexEnter) - lb_N(indexEnter);

            % Find the index that gives the tightest upper bound
            if ~isempty(zind)
                % Harris' test: Pass 1
                bounds = ub_B(zind) + tolprt;
                bounds(indnegz) = lb_B( zind(indnegz) ) - tolprt;
                tB = (bounds - x_B(zind))./ z(zind);
                minUb = min(tB);
                
                % Harris's test: Pass 2
                if minUb <= ub_delta
                    bounds = ub_B(zind);
                    bounds(indnegz) = lb_B( zind(indnegz) );
                    tB = (bounds - x_B(zind))./ z(zind);
                    lvcands = (tB <= minUb);
                    % For an explanation of the code below see the first
                    % "Harris' test: Pass 2" above.
                    [pivotMax_unused,lv] = max(abs(z(zind(lvcands))));
                    indlvcands = find(lvcands,lv);
                    lv = indlvcands(end);
                    minUb = max(0,tB(lv));
                    indexLeave = zind(lv);
                end
            end

            if minUb < ub_delta
                delta = minUb;
                sameBasis = false;
            else
                delta = ub_delta;
                sameBasis = true;
            end

            % Update the solution x_N and x_B
            if delta < inf && delta > tol2
                x_N(indexEnter) = x_N(indexEnter) - delta;
                x_B = x_B + delta * z;
            end
        end

        if ( abs( dualf(indexEnter) ) > tol )
            if isinf(delta)
                exitflag = Unbounded;
                % Set the delta to be 1.0e+16 in unbounded case.
                delta = 1.0e+16;
                if dualf(indexEnter) >= tol
                    x_N(indexEnter) = x_N(indexEnter) - delta;
                    x_B = x_B + delta * z;
                elseif dualf(indexEnter) <= -tol
                    x_N(indexEnter) = x_N(indexEnter) + delta;
                    x_B = x_B - delta * z;
                end
                xsol = [x_B; x_N];
                [tmp, order] = sort([basicVarIdx; nonbasicVarIdx]);
                xsol = xsol(order);
                fval = c'*xsol;
                % Note niters is updated but the last iteration is not executed
                % completely, so niters needs to be adjusted by -1.
                niters = niters - 1;
                % Restore the warning states to their original settings
                warning(warningstate1)
                warning(warningstate2)
                return;
            elseif  delta < tol2
                exitflag = Degenerate;
                dcount = nnz( x_B == lb_B | x_B == ub_B );
                if verbosity >= 5
                    disp( sprintf('******%8d degenerating already.', dcount) );
                    disp( sprintf('indexEnter=%d, \t nonbasicVarIdx(indexEnter) = %d', indexEnter, nonbasicVarIdx(indexEnter) ) );
                    disp( sprintf('indexLeave=%d, \t basicVarIdx(indexLeave) = %d', indexLeave, basicVarIdx(indexLeave) ) );
                    disp( sprintf('    [basicVarIdx      c_B      x_B      lb_B      ub_B] = ') );
                    disp(              full([basicVarIdx      c_B      x_B      lb_B      ub_B]) );
                    disp( sprintf('    [nonbasicVarIdx      c_N      x_N      lb_N      ub_N] = ') );
                    disp(              full([nonbasicVarIdx      c_N      x_N      lb_N      ub_N]) );
                end
            end
        else
            % Restore the warning states to their original settings
            warning(warningstate1)
            warning(warningstate2)
            error('optim:simplexphasetwo:WrongEnteringVar', ...
                'The choice of the entering variable is wrong.');
        end

        % If indexLeave = 0, we are performing a bound flip and Basis stays
        % the same. 
        if ~sameBasis
            % Accept the entering candidate if the corresponding pivot
            % entry is not small.
            if abs(z(indexLeave)) > 1e-8
                acceptCandidate = true;
            else
                % Otherwise, remove this candidate from the candidate list
                % and select another one from the remaining candidates
                % by Devex pricing.
                enteringCandidates(ind_enteringCandidates) = [];
                if ~isempty(enteringCandidates)
                    % Devex pricing 
                    [vl,ind_enteringCandidates] = max(abs(dualf(enteringCandidates).^2)./sqrweights(enteringCandidates));
                    indexEnter = enteringCandidates(ind_enteringCandidates);
                    
                else
                    warning(warningstate1)
                    warning(warningstate2)
                    error('optim:simplexphasetwo:NoEnteringCandidates', ...
                        'There is no other proper choice of the entering variable.');

                end
            end
        else
            % If the same basis, we still accept the candidate to continue
            % with the major iterations
            acceptCandidate = true;
        end % end of "if ~sameBasis"
        
    end % while ~acceptCandidate
    
    % Update Basis, Nonbasis; basicVarIdx, nonbasicVarIdx; and c_B, c_N; ub_B, ub_N; lb_B, lb_N; 
    % x_B(indexLeave), x_N(indexEnter) after each iteration on basis
    % Basis(:,indexLeave) <--> Nonbasis(:, indexEnter);
    % basicVarIdx(indexLeave)   <--> nonbasicVarIdx(indexEnter);
    % c_B(indexLeave)     <--> c_N(indexEnter);
    % ub_B(indexLeave)    <--> ub_N(indexEnter);
    % lb_B(indexLeave)    <--> lb_N(indexEnter);
    % x_B(indexLeave)     <--> x_N(indexEnter);
    
    if ~sameBasis
        Nonbasis(:,indexEnter) = Basis(:,indexLeave);
        Basis(:,indexLeave)    = N_indexEnter;
        
        if ~isempty(eta.cind)
            eta.cind(end+1) = indexLeave;
            eta.matrix(:,end+1) = z;
            eta.pivot(end+1) = z(indexLeave);
        else
            eta.cind(1) = indexLeave; 
            eta.matrix(:,1) = z;  
            eta.pivot(1) = z(indexLeave);
        end
        
        % The number of iterations between refactorizations is
        % determined as 8 by the method described in Chvatal (1983),
        % Linear Programming, p. 111.
        if mod(niters,8) == 0
            RefactorNow = true;
        end
        
        swaptmp           = basicVarIdx(indexLeave);
        basicVarIdx(indexLeave) = nonbasicVarIdx(indexEnter);
        nonbasicVarIdx(indexEnter) = swaptmp;
        
        swaptmp           = c_B(indexLeave); 
        c_B(indexLeave)   = c_N(indexEnter);
        c_N(indexEnter)   = swaptmp; 
        
        swaptmp           = ub_B(indexLeave);
        ub_B(indexLeave)  = ub_N(indexEnter);
        ub_N(indexEnter)  = swaptmp;
        
        swaptmp           = lb_B(indexLeave);
        lb_B(indexLeave)  = lb_N(indexEnter);
        lb_N(indexEnter)  = swaptmp;
        
        swaptmp           = x_B(indexLeave);
        x_B(indexLeave)   = x_N(indexEnter);
        x_N(indexEnter)   = swaptmp;

    end % end of "if ~sameBasis"

    
    xsol = [x_B; x_N];
    fval = c_B' * x_B + c_N' * x_N;
    
    % When maxiter reached, print last iteration and set exitflag = 0.
    if (niters ==  maxiter)
        exitflag = ExcdMaxiter;
        [tmp, order] = sort([basicVarIdx; nonbasicVarIdx]);
        xsol = xsol(order);
        
        % Compute the dual infeasibility when niters == maxiter, 
        % only for display iteration purpose
        if (verbosity >=2) 
            if ~sameBasis
                y = c_B' / Basis;
            end
            dualf = c_N - (y * Nonbasis)';
            enterCandidates = ((dualf >= tol) & (x_N - lb_N > tol2)) | ((dualf <= -tol) & (tol2 < ub_N - x_N));
            disp( sprintf('%8d         %12.6g                %12.6g', niters, full(fval), norm(dualf(enterCandidates))) );
        end
    end
    
    if verbosity >= 5
        % For checking purpose 
        fcheck = max (abs( A(:, basicVarIdx)*xsol(1:m) + A(:,nonbasicVarIdx)*xsol(m+1: n) - b ) );  
        if fcheck > 1.0e-8
            disp( sprintf('  Feasibility is broken in this iteration. %e.', fcheck) );       
        end
        ubcheck = max ( [x_B - ub_B; x_N - ub_N] );
        lbcheck = max ( [lb_B - x_B; lb_N - x_N] );
        if (ubcheck > tol2) || (lbcheck > tol2)
            disp( sprintf('  Boundary condition is violated by %e.', max(ubcheck, lbcheck) ) );
        end
    end
    
end % while (niters < maxiter)

% Restore the warning states to their original settings
warning(warningstate1)
warning(warningstate2)

if (niters == maxiter) 
    exitflag = ExcdMaxiter;
end
%END of SIMPLEXPHASETWO

% ========================================================
% Helper functions for the ETA updates
% ========================================================

function acol = forwardeta(etacind,etapivot,etamatrix,acol)
% This function solves the equation E_1*E_2..*E_k*x = acol
k = length(etacind);
for i = 1:k
    p = etacind(i);
    xptmp = acol(p)./etapivot(i);
    acol = acol - etamatrix(:,i).*xptmp;
    acol(p) = xptmp;
end

function arow = backeta(etacind,etapivot,etamatrix,arow)
% This function solves the equation y*E_1*E_2..*E_k = arow
k = length(etacind);
for i = k:-1:1
    p = etacind(i);
    btmp = arow(p);
    arow(p) = 0;
    btmp = btmp - arow*etamatrix(:,i);
    arow(p) = btmp/etapivot(i);
end
