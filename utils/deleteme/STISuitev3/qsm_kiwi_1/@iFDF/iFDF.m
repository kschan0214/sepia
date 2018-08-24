function  res = iFDF(D2)

res.adjoint = 0;
res.mask = D2;
res = class(res,'iFDF');

