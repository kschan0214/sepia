% Normal Equation of the Forward Calculation used in PDF
%   y = dipole_term(W,D,Mask,xx)
% 
%   output
%   y - a background field 
% 
%   input
%   W - noise covariance matrix
%   D - dipole kernel
%   xx - the background dipoles
%
%   Created by Tian Liu in 2009
%   Last modified by Tian Liu on 2013.07.24

function y = dipole_term(W,D,Mask,xx)

x = zeros(size(D));
x(Mask(:) == 0) = xx(1:end);
x(Mask(:) == 1) = 0;

Ax = real(ifftn(D.* fftn(W.*real(ifftn(D.* fftn(x) )))));
y = Ax( Mask(:) == 0);






