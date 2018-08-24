function res = Wavelet3D(wavScale)
res.adjoint = 0;
[res.af res.sf]= db4;
res.wavScale = wavScale;
res = class(res,'Wavelet3D');
