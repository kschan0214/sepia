% At_dct.m
%

function w = At_dct(y, OMEGA, n)

global numAt;

yu = zeros(n,1);
yu(OMEGA) = y;
w = idct(yu);

numAt = numAt + 1;