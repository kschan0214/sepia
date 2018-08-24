function  res = PhiW(x)

res.adjoint = 0;
res.mask = x;
res = class(res,'PhiW');

