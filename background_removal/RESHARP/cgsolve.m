% cgsolve.m
%
% Solve a symmetric positive definite system Ax = b via conjugate gradients.
%
% Usage: [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose, x0)
%
% A - Either an NxN matrix, or a function handle.
%
% b - N vector
%
% tol - Desired precision.  Algorithm terminates when 
%    norm(Ax-b)/norm(b) < tol .
%
% maxiter - Maximum number of iterations.
%
% verbose - If 0, do not print out progress messages.
%    If and integer greater than 0, print out progress every 'verbose' iters.
%
% x0 - initial solution
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%
% Modified for QSM by Ildar Khalidov, WCMC
% Last Modified by Tian Liu on 2011.02.02

function [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose, x0)

matrix_size=size(b);
b=b(:);


if (nargin < 6)
    x = zeros(length(b),1);
else 
    x=x0(:);
end

if (nargin < 5), verbose = 1; end

implicit = isa(A,'function_handle');

if (nargin < 6)
    r = b;
else
    if (implicit), r = b - A(reshape(x,matrix_size)); r=r(:);  else, r = b - A*x;  end
end
d = r;
delta = r'*r;
delta0 = b'*b;
numiter = 0;
bestx = x;
bestres = sqrt(delta/delta0); 
while ((numiter < maxiter) & (delta > tol^2*delta0))

  % q = A*d
  if (implicit), q = A(reshape(d,matrix_size)); q=q(:);  else, q = A*d;  end
 
  alpha = delta/(d'*q);
  x = x + alpha*d;
  
  if (mod(numiter+1,50) == 0)
    % r = b - Aux*x
    if (implicit), r = b - reshape(A(reshape(x,matrix_size)),size(b));  else, r = b - A*x;  end
  else
    r = r - alpha*q;
  end
  
  deltaold = delta;
  delta = r'*r;
  beta = delta/deltaold;
  d = r + beta*d;
  numiter = numiter + 1;
  if (sqrt(delta/delta0) < bestres)
    bestx = x;
    bestres = sqrt(delta/delta0);
  end    
  
  if ((verbose) & (mod(numiter,verbose)==0))
    disp(sprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e', ...
      numiter, bestres, sqrt(delta/delta0)));
  end
  
end

if (verbose)
  disp(sprintf('cg: Iterations = %d, best residual = %14.8e', numiter, bestres));
end
% x = reshape(bestx,matrix_size);
x = reshape(x,matrix_size);
res = bestres;
iter = numiter;

