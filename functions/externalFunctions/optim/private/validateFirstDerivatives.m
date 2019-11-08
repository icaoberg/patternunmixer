function validateFirstDerivatives(funfcn,confcn,x,fval,cIneq,cEq,grad,JacCineqTrans,JacCeqTrans, ...
                                          lb,ub,fscale,options,sizes,varargin)
% validateFirstDerivatives Helper function that validates first derivatives of
% objective, nonlinear inequality, and nonlinear equality gradients against
% finite differences. The finite-difference calculation is done according to
% options.FinDiffType.

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2007/12/10 21:50:57 $

tol = 1e-6; % component-wise relative difference in gradients checked against this tolerance
mNonlinIneq = sizes.mNonlinIneq;  
mNonlinEq = sizes.mNonlinEq;

if strcmpi(options.GradObj,'on') 
    % Finite-difference gradient already scaled appropriately since computed with scaled quantities
    grad_fd = finDiffObj(funfcn{3},x,fval,lb,ub,fscale,options,sizes,varargin{:});
    
    % Unscale gradients
    if strcmpi(options.ScaleProblem,'obj-and-constr')
        grad = grad / fscale.obj;        
        grad_fd = grad_fd / fscale.obj;
    end

    % Vector of objective gradient relative error
    gradError = abs(grad_fd - grad);    
    [maxDiff,grad_idx] = max(gradError);
    fprintf('Objective function derivatives:\n')
    fprintf('Maximum discrepancy between derivatives: %g\n',maxDiff)
    if any(gradError./max(1.0,grad) > tol)
        fprintf('Warning: user-supplied and %s finite-difference derivatives do\n',options.FinDiffType)
        fprintf(' not match within %g relative tolerance.\n',tol)
        fprintf('Maximum difference occurs in element %i of gradient:\n',grad_idx)
        fprintf('  User-supplied gradient:     % g\n',grad(grad_idx))
        fprintf('  Finite-difference gradient: % g\n',grad_fd(grad_idx))
        disp('Strike any key to continue or Ctrl-C to abort.')
        pause
        disp(' ') % blank line
    end
end

% If there are nonlinear constraints and their derivatives are provided,
% validate them
if strcmpi(options.GradConstr,'on') && sizes.mNonlinIneq + sizes.mNonlinEq > 0
    % Finite-diffence Jacobians already scaled appropriately since computed with scaled quantities
    [JacCineqTrans_fd,JacCeqTrans_fd] = ...
        finDiffConstr(confcn{3},x,-cIneq(:),cEq(:),lb,ub,fscale,options,sizes,varargin{:});

    if sizes.mNonlinIneq > 0
        % Unscale Jacobian
        if strcmpi(options.ScaleProblem,'obj-and-constr')
            JacCineqTrans = JacCineqTrans * spdiags(1.0./fscale.cIneq,0,mNonlinIneq,mNonlinIneq);
            JacCineqTrans_fd = JacCineqTrans_fd * spdiags(1.0./fscale.cIneq,0,mNonlinIneq,mNonlinIneq);
        end

        % Matrix of nonlinear inequality constraint gradient relative error
        % JacCineqTrans_fd is full so JacCineqError will be full - store it as a full matrix        
        JacCineqError = full(abs(JacCineqTrans - JacCineqTrans_fd));
        [maxDiff,i,j] = findRowColIndicesOfMaxElement(JacCineqError);
        fprintf('Nonlinear inequality constraint derivatives:\n')
        fprintf('Maximum discrepancy between derivatives: %g\n',maxDiff)
        if any(any( JacCineqError./max(1.0,JacCineqTrans) > tol ))
            fprintf('Warning: user-supplied and %s finite-difference derivatives do\n',options.FinDiffType)
            fprintf(' not match within %g relative tolerance.\n',tol)
            fprintf('Maximum difference occurs in element (%i,%i):\n',i,j)
            fprintf('  User-supplied constraint gradient:     % g\n',JacCineqTrans(i,j))
            fprintf('  Finite-difference constraint gradient: % g\n',JacCineqTrans_fd(i,j))
            disp('Strike any key to continue or Ctrl-C to abort.')
            pause
            disp(' ') % blank line
        end
    end
    
    if sizes.mNonlinEq > 0
        % Unscale Jacobian
        if strcmpi(options.ScaleProblem,'obj-and-constr')
            JacCeqTrans = JacCeqTrans * spdiags(1.0./fscale.cEq,0,mNonlinEq,mNonlinEq);
            JacCeqTrans_fd = JacCeqTrans_fd * spdiags(1.0./fscale.cEq,0,mNonlinEq,mNonlinEq);
        end
        % Matrix of nonlinear equality constraint gradient relative error
        % JacCeqTrans_fd is full so JacCeqError will be full - store it as a full matrix
        JacCeqError = full(abs(JacCeqTrans - JacCeqTrans_fd)); 
        [maxDiff,i,j] = findRowColIndicesOfMaxElement(JacCeqError);
        fprintf('Nonlinear equality constraint derivatives:\n')
        fprintf('Maximum discrepancy between derivatives: %g\n',maxDiff)
        if any(any( JacCeqError./max(1.0,JacCeqTrans) > tol ))
            fprintf('Warning: user-supplied and %s finite-difference derivatives do\n',options.FinDiffType)
            fprintf(' not match within %g relative tolerance.\n',tol)
            fprintf('Maximum difference occurs in element (%i,%i):\n',i,j)
            fprintf('  User-supplied constraint gradient:     % g\n',JacCeqTrans(i,j))
            fprintf('  Finite-difference constraint gradient: % g\n',JacCeqTrans_fd(i,j))
            fprintf('Strike any key to continue or Ctrl-C to abort.')
            pause
            disp(' ') % blank line
        end
    end
end

%-------------------------------------------------------------------------
function [maxVal,i,j] = findRowColIndicesOfMaxElement(A)
% Helper function that finds indices (i,j) of the maximum element of matrix A.
% It also returns the maximum element, maxVal

[col_max_val,row_idx] = max(A);
[maxVal,col_idx] = max(col_max_val);
i = row_idx(col_idx);
j = col_idx;





