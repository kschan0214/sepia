% A_dct.m
%

function y = A_dct(x, OMEGA)

global numA;

yu = dct(x);
y = yu(OMEGA);

numA = numA + 1;