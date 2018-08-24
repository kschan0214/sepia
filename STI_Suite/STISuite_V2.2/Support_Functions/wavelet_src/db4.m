function [af, sf] = db4

% Farras nearly symmetric filters for orthogonal
% 2-channel perfect reconstruction filter bank
%
% USAGE:
%    [af, sf] = farras
% OUTPUT:
%    af - analysis filters
%    sf - synthesis filters
% REFERENCE:
%    A. F. Abdelnour and I. W. Selesnick. 
%    "Nearly symmetric orthogonal wavelet bases",
%    Proc. IEEE Int. Conf. Acoust., Speech,
%    Signal Processing (ICASSP), May 2001.
% See afb, dwt.
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
%last change by Shuai 2010.6.29 for create db4
af = [
-0.010597401784997278        -0.23037781330885523 
0.032883011666982945         0.71484657055254153  
0.030841381835986965         -0.63088076792959036 
-0.18703481171888114         -0.027983769416983849
-0.027983769416983849        0.18703481171888114  
0.63088076792959036          0.030841381835986965 
0.71484657055254153          -0.032883011666982945
0.23037781330885523          -0.010597401784997278
];
 
sf = af(end:-1:1, :);