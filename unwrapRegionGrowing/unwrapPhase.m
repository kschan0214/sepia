% Unwrap phase using a region growing algorithm
%   [iFreq] = unwrapPhase(iMag, iFreq_raw, matrix_size)
% 
%   output
%   iFreq - unwrapped phase
%
%   input
%   iMag - quality map reflecing reliability of each point
%   iFreq_raw - A wrapped field map
%   matrix_size - dimension of the field of view
%   
%   C version created by Tian Liu in 2010.
%   MATLAB Wrapped created by Dong Zhou(zhou.dong@gmail.com) on 2013.06.12
%   Modified by Dong Zhou on 2013.06.24
%   Last modified by Tian Liu on 2013.07.23


function [iFreq] = unwrapPhase(iMag, iFreq_raw, matrix_size)
    matrix_size  = double(matrix_size);
    iMag = double(iMag);
    iFreq_raw = double(iFreq_raw);
    iFreq = mexUnwrap(iMag(:),iFreq_raw(:),matrix_size(1),matrix_size(2),matrix_size(3));
    iFreq = reshape(iFreq,matrix_size);
end
