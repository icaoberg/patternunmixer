% Optimization Toolbox
% Version 4.0 (R2008a) 23-Jan-2008 
%
% Nonlinear minimization of functions.
%   fminbnd      - Scalar bounded nonlinear function minimization.
%   fmincon      - Multidimensional constrained nonlinear minimization.
%   fminsearch   - Multidimensional unconstrained nonlinear minimization, 
%                  by Nelder-Mead direct search method.
%   fminunc      - Multidimensional unconstrained nonlinear minimization.
%   fseminf      - Multidimensional constrained minimization, semi-infinite 
%                  constraints.
%
% Nonlinear minimization of multi-objective functions.
%   fgoalattain  - Multidimensional goal attainment optimization 
%   fminimax     - Multidimensional minimax optimization.
%        
% Linear least squares (of matrix problems).
%   lsqlin       - Linear least squares with linear constraints.
%   lsqnonneg    - Linear least squares with nonnegativity constraints.
%
% Nonlinear least squares (of functions).
%   lsqcurvefit  - Nonlinear curvefitting via least squares (with bounds).
%   lsqnonlin    - Nonlinear least squares with upper and lower bounds.
%
% Nonlinear zero finding (equation solving).
%   fzero        - Scalar nonlinear zero finding.
%   fsolve       - Nonlinear system of equations solve (function solve).
%
% Minimization of matrix problems.
%   bintprog     - Binary integer (linear) programming.
%   linprog      - Linear programming.
%   quadprog     - Quadratic programming.
%
% Controlling defaults and options.
%   optimset     - Create or alter optimization OPTIONS structure. 
%   optimget     - Get optimization parameters from OPTIONS structure. 
%
%
% Graphical user interface and plot routines
%   optimtool                   - Optimization Toolbox Graphical User 
%                                 Interface
%   optimplotconstrviolation    - Plot max. constraint violation at each 
%                                 iteration
%   optimplotfirstorderopt      - Plot first-order optimality at each 
%                                 iteration
%   optimplotresnorm            - Plot value of the norm of residuals at
%                                 each iteration
%   optimplotstepsize           - Plot step size at each iteration

% Internally Used Utility Routines
%
%   Large-scale preconditioners
%   aprecon      - default banded preconditioner for least-squares problems
%   hprecon      - default banded preconditioner for minimization problems
%   pceye        - diagonal preconditioner based on scaling matrices
%
%   Large-scale utility routines
%   fzmult       - Multiplication with fundamental nullspace basis
%   gangstr      - Zero out 'small' entries subject to structural rank
%   sfminbx      - Nonlinear minimization with box constraints
%   sfminle      - Nonlinear minimization with linear equalities
%
%   Cubic interpolation routines
%   cubici1      - interpolates 2 pts and gradients to estimate minimum
%   cubici2      - interpolates 3 points and 1 gradient
%   cubici3      - interpolates 2 pts and gradients to find step and min
%
%   Quadratic interpolation routines
%   quadi        - interpolates three points to estimate minimum
%
%   Semi-infinite utility routines
%   semicon      - translates semi-infinite constraints to constrained problem
%   findmax      - interpolates the maxima in a vector of data 
%   v2sort       - sorts two vectors and then removes missing elements
%
%   Goal-attainment utility routines
%   goalfun      - translates goal-attainment problem to constrained problem
%   goalcon      - translates constraints in goal-attainment problem
%
%   Graphical user interface utility routines
%   optimguiswitchyard  - Switchyard routine for optimtool
%   optimtooloutput     - OutputFcn used to interact between the GUI and 
%                         solvers
%
%   Other
%   color                 - column partition for sparse finite differences
%   graderr               - used to check gradient discrepancy 
%   searchq               - line search routine for lsqnonlin, lsqcurvefit,
%                           and fsolve functions
%   sfd                   - sparse Hessian via finite gradient differences
%   sfdnls                - sparse Jacobian via finite differences
%   optimoptions          - display option names and values for the options
%                           in the Optimization Toolbox but not in MATLAB
%   optimoptioncheckfield - check validity of field contents
%   optimoptiongetfields  - fieldnames of options in Optimization Toolbox 
%                           not in MATLAB
%   functiontostring      - converts the function name to a string for 
%                           display

%   Copyright 1990-2008 The MathWorks, Inc.
%   Generated from Contents.m_template revision 1.1.6.4  $Date: 2007/12/10 21:49:44 $
