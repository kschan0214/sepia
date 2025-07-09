%% [RDF] = Wrapper_BFR_PDF(totalField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)
%
% Input
% --------------
% totalField    : total field map (background + tissue fields), in Hz
% mask          : signal mask
% matrixSize    : size of the input image
% voxelSize     : spatial resolution of each dimension of the data, in mm
% algorParam    : structure contains fields with algorithm-specific parameter(s)
% headerAndExtraData : structure contains extra header info/data for the algorithm
%
% Output
% --------------
% RDF           : local field map
%
% Description: This is a wrapper function to access PDF for SEPIA
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 8 March 2020 (v0.8.0)
% Date last modified:
%
%
function [unwrappedField] = Wrapper_Unwrap_None(wrappedField,mask,matrixSize,voxelSize,algorParam, headerAndExtraData)

unwrappedField = wrappedField;

end