%% function [Ex,Ey,Ez,EtE] = GradientOperatorKspace(matrixSize)
%
% Input
% --------------
%   matrixSize      : matrix size of chi
%
% Output
% --------------
%   Ex              : gradient in x-direction
%   Ey              : gradient in y-direction
%   Ez              : gradient in z-direction
%   EtE             : magnitude of gradient 
%
% Description: Gradient operator in K-space
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 15 July 2017
% Date last modified:
%
function [Ex,Ey,Ez,EtE] = GradientOperatorKspace(matrixSize)
[k1,k2,k3] = ndgrid(-matrixSize(1)/2:matrixSize(1)/2-1, ...
                    -matrixSize(2)/2:matrixSize(2)/2-1, ...
                    -matrixSize(3)/2:matrixSize(3)/2-1);
% KC: gradient terms in fourier space
Ex = fftshift(1 - exp(2i*pi .* k1 / matrixSize(1)));
Ey = fftshift(1 - exp(2i*pi .* k2 / matrixSize(2)));
Ez = fftshift(1 - exp(2i*pi .* k3 / matrixSize(3)));
EtE = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;
end