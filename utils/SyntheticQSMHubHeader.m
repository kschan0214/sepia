%% [b0,b0dir,voxelSize,matrixSize,te,dte,CF]=SyntheticQSMHubHeader(multiecho,varargin)
%
% Usage:
%
% Input
% --------------
%
% Output
% --------------
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 11 April 2018
% Date last modified:
%
%
function [b0,b0dir,voxelSize,matrixSize,te,dte,CF]=SyntheticQSMHubHeader(multiecho,varargin)
gyro = 42.57747892;

[matrixSize(1),matrixSize(2),matrixSize(3),nt] = size(multiecho);

[b0,b0dir,voxelSize,te] = process_options(varargin,'b0',3,'b0dir',[0,0,1],'voxel',[1,1,1],'te',linspace(1e-3,30e-3,nt));

dte = te(2)-te(1); %s

CF = gyro * 1e6 * b0;
            
end