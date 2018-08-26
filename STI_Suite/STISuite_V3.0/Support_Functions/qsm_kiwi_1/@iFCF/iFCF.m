function  res = iFCF(C)

res.adjoint = 0;
res.mask = C;
res = class(res,'iFCF');

