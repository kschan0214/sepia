%% function [resnorm,residual] = ResidualGivenQSMLocalFieldDipolekernel(chi,localField,kernel)
%
% Description: compute the residual and residual norm given the local
%              field, dipole kernel and QSM
%
% Input
% -----
%   chi             : QSM
%   localField      : local field perturbatios
%   kernel          : dipole kernel
%
% Output
% ______
%   resnorm     	: residual norm
%   residual        : residual map
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 26 March 2017
% Date last modified: 28 June 2017
%
function [resnorm,residual] = ResidualGivenQSMLocalFieldDipolekernel(chi,localField,kernel)

residual = ifftn(kernel.*fftn(chi)) - localField;
resnorm = norm(residual(:),2);

end