function res = mtimes(a,b)

res = ifftnc(a.mask.*fftnc(b));

