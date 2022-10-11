%% function u = ARLO(y,x)
%
% Description: Fast monoexponential fitting by auto-regression
% Samples have to be evenly sampled (i.e. even spacing)
% At least 3 samples are needed
% Assuming time series in the last dimension
% ref: Pei et al. MRM 73:843-850(2015)
% e.g. for function y=exp(-x/u), u can be estimated by ARLO
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 October, 2016
% Date last modified:
%
function u = ARLO(y,x)
% ensure y is real
% y=abs(y);

% get dimension of y
ndim = ndims(y);
if ndim==2 && size(y,2)==1
    y = permute(y,[2 1]);
end
matrixSize = size(y);
N = matrixSize(end);

% reshape y s.t. new y = [all y, time]
y = reshape(y,[numel(y)/matrixSize(end)],N);

% get the spacing
dx = x(2)-x(1);

% get sum of signal for i and delta i
Si = zeros(size(y,1),N-2);
deltai = zeros(size(y,1),N-2);
for k=1:N-2
    [Si(:,k), deltai(:,k)]= Simpson(y(:,k:k+2),dx);
end

% analytical solution for minimiser to obatain u
a = sum(Si.^2,2);
b = sum(Si.*deltai,2);
u = (a + (dx/3)*b)./((dx/3)*a + b);

% reshape u based on input dimension
u = reshape(u,[matrixSize(1:end-1)]);

end

% Quadratic approximation of Simpson rule's in 4th order accuracy when J=2
function [Si, deltai] = Simpson(y,dx)
Si = (dx/3) * (y(:,1) + 4*y(:,2) + y(:,3));
deltai = y(:,1)-y(:,3);
end