function res = mtimes(a,b)
if isa(a,'myWavelet') == 0
    error('In  A.*B only A can be Wavelet operator');
end

if a.adjoint
    res = idwt3D(a.wavScale,a.sf,real(b)) + 1i*idwt3D(a.wavScale,a.sf,imag(b));
else
    res = dwt3D(real(b),a.wavScale,a.af) + 1i*dwt3D(imag(b),a.wavScale,a.af);
end


