function NewImag=applyxfm(imag,XFM)

rot=XFM(1:3,1:3);
rot=rot/det(rot)^(1/3);
rot(4,4)=1;
NewImag = affine(imag(end:-1:1,:,:), rot);
NewImag = NewImag(end:-1:1,:,:);