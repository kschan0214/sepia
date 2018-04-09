%  SpaRSA version 2.0, December 31, 2007
% modified by H.Wei Duke 2014
%  
%  This function solves the convex problem 
% 
%  arg min_x = 0.5*|| y - A x ||_2^2 + tau phi(x)
% 
%   d'*A'*(A x - y) + 0.5*alpha*d'*d 
% 
%  where alpha is obtained from a BB formula. In a monotone variant, alpha is
%  increased until we see a decreasein the original objective function over
%  this step. 
% 
%  -----------------------------------------------------------------------
%  Copyright (2007): Mario Figueiredo, Robert Nowak, Stephen Wright
%  
%  ----------------------------------------------------------------------
%  
% 
%   ===== Required inputs =============
% 
%   y: 1D vector or 2D array (image) of observations
%      
%   A: if y and x are both 1D vectors, A can be a 
%      k*n (where k is the size of y and n the size of x)
%      matrix or a handle to a function that computes
%      products of the form A*v, for some vector v.
%      In any other case (if y and/or x are 2D arrays), 
%      A has to be passed as a handle to a function which computes 
%      products of the form A*x; another handle to a function 
%      AT which computes products of the form A'*x is also required 
%      in this case. The size of x is determined as the size
%      of the result of applying AT.
% 
%   tau: regularization parameter (scalar)
% 
%   ===== Optional inputs =============
% 
%   
%   'AT'    = function handle for the function that implements
%             the multiplication by the conjugate of A, when A
%             is a function handle. If A is an array, AT is ignored.
% 
%   'Psi'   = handle to the denoising function, that is, to a function
%             that computes the solution of the densoing probelm 
%             corresponding to the desired regularizer. That is, 
%             Psi(y,tau) = arg min_x (1/2)*(x - y)^2 + tau phi(x).
%             Default: in the absence of any Phi given by the user,
%             it is assumed that phi(x) = ||x||_1 thus 
%             Psi(y,tau) = soft(y,tau)
%             Important: if Psi is given, phi must also be given,
%                        so that the algorithm may also compute
%                        the objective function.
% 
%   'StopCriterion' = type of stopping criterion to use
%                     0 = algorithm stops when the relative 
%                         change in the number of non-zero 
%                         components of the estimate falls 
%                         below 'ToleranceA'
%                     1 = stop when the relative 
%                        change in the objective function 
%                        falls below 'ToleranceA'
%                     2 = stop when relative duality gap 
%                        falls below 'ToleranceA'
%                     3 = stop when LCP estimate of relative
%                        distance to solution falls below ToleranceA
%                     4 = stop when the objective function 
%                         becomes equal or less than toleranceA.
%                     5 = stop when the norm of the difference between 
%                         two consecutive estimates, divided by the norm
%                         of one of them falls below toleranceA
%                     Default = 2
% 
%   'ToleranceA' = stopping threshold; Default = 0.01
%  
%   'Debias'     = debiasing option: 1 = yes, 0 = no.
%                  Default = 0.
% 
%   'ToleranceD' = stopping threshold for the debiasing phase:
%                  Default = 0.0001.
%                  If no debiasing takes place, this parameter,
%                  if present, is ignored.
% 
%   'MaxiterA' = maximum number of iterations allowed in the
%                main phase of the algorithm.
%                Default = 1000
% 
%   'MiniterA' = minimum number of iterations performed in the
%                main phase of the algorithm.
%                Default = 5
% 
%   'MaxiterD' = maximum number of iterations allowed in the
%                debising phase of the algorithm.
%                Default = 200
% 
%   'MiniterD' = minimum number of iterations to perform in the
%                debiasing phase of the algorithm.
%                Default = 5
% 
%   'Initialization' must be one of {0,1,2,array}
%                0 -> Initialization at zero. 
%                1 -> Random initialization.
%                2 -> initialization with A'*y.
%            array -> initialization provided by the user.
%                Default = 0;
% 
%   'BB_variant' specifies which variant of Barzila-Borwein to use, or not.
%                0 -> don't use a BB rule - instead pick the starting alpha
%                based on the successful value at the previous iteration
%                1 -> standard BB choice  s'r/s's
%                2 -> inverse BB variant r'r/r's
%                Default = 1
% 
%   'BB_cycle' specifies the cycle length  - the number of iterations between
%              recalculation of alpha. Requires integer value at least
%              1. Relevant only if a **nonmonotone BB rule** is used 
%              (BB_variant = 1 or 2 and Monotone=0).
%              Default = 1
% 
%   'Monotone' =  enforce monotonic decrease in f, or not? 
%                any nonzero -> enforce monotonicity (overrides 'Safeguard')
%                0 -> don't enforce monotonicity.
%                Default = 0;
% 
%   'Safeguard' = enforce a "sufficient decrease" over the largest
%                objective value of the past M iterations.
%                any nonzero -> safeguard
%                0 -> don't safeguard
%                Default = 0.
% 
%   'M'        = number of steps to look back in the safeguarding process.
%                Ignored if Safeguard=0 or if Monotone is nonzero.
%                (positive integer. Default = 5)
% 
%   'sigma'    = sigma value used in Safeguarding test for sufficient 
%                decrease. Ignored unless 'Safeguard' is nonzero. Must be
%                in (0,1). Drfault: .01.
% 
%   'Eta'      = factor by which alpha is multiplied within an iteration,
%                until a decrease in the objective function is
%                obtained.
%                Default = 2;
% 
%   'Alpha_factor' = factor by which to reduce the successful value of
%                 alpha at iteration k, to give the first value of alpha
%                 to be tried at iteration k+1.
%                 If a Barzilai-Borwein rule is specified (BB_variant > 0), 
%                 this parameter is ignored.
%                 Default = 0.8;
% 
%   'Continuation' = Continuation or not (1 or 0) 
%                    Specifies the choice for a continuation scheme,
%                    in which we start with a large value of tau, and
%                    then decrease tau until the desired value is 
%                    reached. At each value, the solution obtained
%                    with the previous values is used as initialization.
%                    Default = 0
% 
%  'ContinuationSteps' = Number of steps in the continuation procedure;
%                        ignored if 'Continuation' equals zero.
%                        If -1, an adaptive continuation procedure is used.
%                        Default = -1.
%  
%  'FirstTauFactor'  = Initial tau value, if using continuation, is
%                      obtained by multiplying the given tau by 
%                      this factor. This parameter is ignored if 
%                      'Continuation' equals zero or 
%                      'ContinuationSteps' equals -1.
%                      Default = 10.
% 
%   'True_x' = if the true underlying x is passed in 
%                 this argument, MSE plots are generated.
% 
%   'AlphaMin' = the alphamin parameter of the BB method.
%                Default = 1e-30;
% 
%   'AlphaMax' = the alphamax parameter of the BB method.
%                Default = 1e30;
% 
%   'Verbose'  = work silently (0) or verbosely (1)
% 
%  ===================================================  
%  ============ Outputs ==============================
%    x = solution of the main algorithm
% 
%    x_debias = solution after the debiasing phase;
%                   if no debiasing phase took place, this
%                   variable is empty, x_debias = [].
% 
%    objective = sequence of values of the objective function
% 
%    times = CPU time after each iteration
% 
%    debias_start = iteration number at which the debiasing 
%                   phase started. If no debiasing took place,
%                   this variable is returned as zero.
% 
%    mses = sequence of MSE values, with respect to True_x,
%           if it was given; if it was not given, mses is empty,
%           mses = [].
%  ========================================================
%
