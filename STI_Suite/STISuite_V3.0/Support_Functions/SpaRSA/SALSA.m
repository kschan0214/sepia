function [z, objective, distance, times, mses] = ...
         SALSA(y,OMEGA,tau,n,varargin)

%%%% Parameter values
%% CONVOLUTIONFILTER    isAconvfilter = 1 => A is the convolution filter, 
%%                      of the same size as the image and centered
%%                      isAconvfilter = 0 => A is either a matrix of
%%                      function handle. (default: 0)
%%
%% WAVELET              iswaveletrepr = 1,2 => the problem is solved for the
%%                      wavelet representation of the image. 2 ->
%%                      redundant, 1 -> orthogonal. (default: 0)

%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 4
     error('Wrong number of required parameters');
end

%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
stopCriterion = 1;
compute_mse = 0;
tolA = 0.01;
maxiter = 10000; % outer loop iterations
init = 0;
AT = 0;
mu = 1e-7;
inneriters = 1;
psi_ok = 0;
tolA = 0.001;
isAconvfilter = 0;
iswaveletrepr = 0;
isAmatrix = 0;
verbose = 1;
definedbasis = 0;
definedbasistranspose = 0;
mses = [];

%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch upper(varargin{i})
     case 'BASIS'
         definedbasis = 1;
       W = varargin{i+1};   
     case 'BASISTRANSPOSE'
         definedbasistranspose = 1;
       WT = varargin{i+1};   
     case 'CONVOLUTIONFILTER'
       isAconvfilter = varargin{i+1};
     case 'WAVELET'
       iswaveletrepr = varargin{i+1};
     case 'PSI'
       psi_function = varargin{i+1};
     case 'PHI'
       phi_function = varargin{i+1};
     case 'MU'
       mu = varargin{i+1};
     case 'STOPCRITERION'
       stopCriterion = varargin{i+1};
     case 'TOLERANCEA'       
       tolA = varargin{i+1};
     case 'INNERITERS'
         inneriters = varargin{i+1};
     case 'MAXITERA'
       maxiter = varargin{i+1};
     case 'INITIALIZATION'
       if prod(size(varargin{i+1})) > 1   % we have an initial x
	 init = 33333;    % some flag to be used below
	 x = varargin{i+1};
       else 
	 init = varargin{i+1};
       end
     case 'TRUE_X'
       compute_mse = 1;
       true = varargin{i+1};
       if prod(double((size(true) == size(y))))
           plot_ISNR = 1;
       end
     case 'AT'
       AT = varargin{i+1};
     case 'VERBOSE'
         verbose = varargin{i+1};
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized option: ''' varargin{i} '''']);
    end;
  end;
end
%%%%%%%%%%%%%%

if (sum(stopCriterion == [0 1 2 3])==0)
   error(['Unknown stopping criterion']);
end


% Precompute A'*y since it'll be used a lot
ATy = At_dct(y, OMEGA, n);


% set functions
psi_function = @(x,tau) soft(x,tau);

phi_function = @(x) sum(abs(x(:))); 
phi_l1 = 1;


% initializing
% initialize at zero
x = 0*ATy;
z = 0*ATy;
b = z;
xprev = x;
zprev = realmax*ones(size(z));
muinv = 1./mu;
threshold = tau./mu;
criterion(1) = 1;

% Compute and store initial value of the objective function
resid =  y - A_dct(x,OMEGA);
prev_f = 0.5*(resid(:)'*resid(:)) + tau*phi_function(x);

if verbose
    fprintf('Initial value of objective function = %3.3g\n',prev_f)
end

% start the clock
t0 = cputime;

times(1) = 0;
objective(1) = prev_f;

% Wx = W(x);
% if compute_mse
%     Wxtrue = W(true);
%    %mses(1) = norm(Wx(:)-Wxtrue(:),2)^2 /numel(Wx);
%    mses(1) = sum(sum((x-true).^2))/numel(x);
% end

mu = tau;

for outer = 1:maxiter
    
    %if (norm(z(:)-x(:),2) >= norm(zprev(:)-xprev(:),2))&&(outer > 5)
    mu = min(mu * 1.15,1);
    
    %end
  
    xprev = x;
    zprev = z;
    
    muinv = 1/mu;
    threshold = tau/mu;
    mu_factor = 1/(1+mu); 
    
   for inner = 1:1
        
           xp = x;
           z = psi_function( x - b, threshold);
           r = ATy + mu*(z+b);
           x = muinv*( r - mu_factor*At_dct( A_dct(r,OMEGA) ,OMEGA,n));
                   
   end

    b = b + (z - x);
    
   resid = y - A_dct(z,OMEGA);
   objective(outer+1) = 0.5*(resid(:)'*resid(:)) + tau*phi_function(z);
    
    if compute_mse
        Wx = W(x);
        %mses(outer+1) = norm(Wx(:)-Wxtrue(:),2)^2 /numel(Wx);
        mses(outer+1) =  sum(sum((x-true).^2))/numel(x);
    end
    
    distance(outer) = norm(x(:)-z(:),2);
                      %sqrt(norm(x(:),2)^2 + norm(z(:),2)^2);

    
    if (outer>1)
        % take no less than miniter and no more than maxiter iterations
        switch stopCriterion
            case 1
                % compute the stopping criterion based on the relative
                % variation of the objective function.
                %criterion(outer) = abs(obj(outer)-obj(outer-1))/obj(outer-1);
                criterion(outer) = abs(objective(outer+1)-objective(outer))/objective(outer);
            case 2
                % compute the stopping criterion based on the relative
                % variation of the estimate.
                criterion(outer) = abs(norm(x(:)-xprev(:))/norm(x(:)));
            case 3
                % continue if not yet reached target value tolA
                criterion(outer) = objective(outer+1);
            otherwise
                error(['Unknown stopping criterion']);
        end
        
         if ( criterion(outer) <= tolA )
             if verbose
                 fprintf('iter = %d, obj = %3.3g, ||x-z|| = %g, stop_cri = %3.3g, (goal = %3.3g ), mu = %g\n', outer, objective(outer+1), distance(outer), criterion(outer), tolA,mu)
                 fprintf('Convergence reached for Bregman iterative procedure.\n')
             end
             times(outer+1) = cputime - t0;
             break;
         end
    end
    
    if verbose
        %fprintf('iter = %d, obj = %3.3g, ||x-z|| = %3.3g, stop criterion = %3.3g\n', outer, obj(outer), distance(outer), criterion(outer))
        %fprintf('iter = %d, obj = %3.3g, ||x-z|| = %3.3g, stop criterion = %3.3g\n', outer, objective(outer+1), distance(outer), criterion(outer))
       fprintf('iter = %d, obj = %3.3g, ||x-z|| = %g, stop_cri = %3.3g, (goal = %3.3g ), mu = %g\n', outer, objective(outer+1), distance(outer), criterion(outer), tolA,mu)
    end
    
    times(outer+1) = cputime - t0;
end

times(outer) = cputime - t0;

tfinal = cputime - t0;

if verbose
    fprintf('Breg CPU time = %3.3g seconds, iters = %d\n', tfinal, outer)

    if compute_mse
        mse = norm(Wxtrue(:)-Wx(:),2)^2 /numel(Wx);
        fprintf('MSE = %3.3g\n', mse)
        ISNR = 10*log10(  norm(y(:)-Wxtrue(:),2)^2 / (mse*numel(Wx)) );
        fprintf('MSE = %3.3g, ISNR = %3.3g dB\n', mse, ISNR)
    end
end
