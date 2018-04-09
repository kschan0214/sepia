function  res = suscepMap(mask,D2)

res.mask = mask;
res.adjoint = 0;
res.D2 = D2;
res = class(res,'suscepMap');