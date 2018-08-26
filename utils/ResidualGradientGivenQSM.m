%% function [resnorm, residual] = ResidualGradientGivenQSM(chi,matrixSize)
%
% Description: compute residual adn residual norm of gradient given QSM
%
% Input
% -----
%   chi             : QSM
%   matrixSize      : image matrix size
%
% Output
% ______
%   resnorm       	: residual norm
%   residual     	: residual map
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 26 March 2017
% Date last modified: 28 June 2017
%
function [resnorm, residual] = ResidualGradientGivenQSM(chi,matrixSize)
% defining gradient operators in k-space
[k1,k2,k3] = ndgrid(-matrixSize(1)/2:matrixSize(1)/2-1, ...
                    -matrixSize(2)/2:matrixSize(2)/2-1, ...
                    -matrixSize(3)/2:matrixSize(3)/2-1);
% KC: gradient terms in fourier space
Ex = fftshift(1 - exp(2i*pi .* k1 / matrixSize(1)));
Ey = fftshift(1 - exp(2i*pi .* k2 / matrixSize(2)));
Ez = fftshift(1 - exp(2i*pi .* k3 / matrixSize(3)));

% KC: Gradient norm
residual(:,:,:,1) = ifftn(Ex.*fftn(chi)) ;
residual(:,:,:,2) = ifftn(Ey.*fftn(chi)) ;
residual(:,:,:,3) = ifftn(Ez.*fftn(chi)) ;
resnorm = norm(residual(:),2);
end