function res = mtimes(a,b)

if a.adjoint   
    %res = ifftnc(a.D2.*fftnc(b)).*a.mask;
   % res = ifftnc(a.D2.*fftnc(b.*a.mask)).*a.mask;
   %  res = ifftnc(a.D2.*fftnc(b));  
   % res = (a.D2.*fftnc(b)); % to work in kspace
 %  res = ifftnc((a.D2.^2).*fftnc(b));
  % res = ifftnc((a.D2.^2).*fftnc(b.*a.mask)).*a.mask;
 % res = ifftnc((a.D2.^2).*fftnc(b.*a.mask));
  res = ifftnc(b.*(a.D2)).*a.mask;
  %res = ifftnc(b).*a.mask;
else
    %res = ifftnc(a.D2.*fftnc(b.*a.mask));
   % res = ifftnc(a.D2.*fftnc(b.*a.mask)).*a.mask;
  % res = ifftnc(a.D2.*(b)); % to work in kspace
  % res = ifftnc((a.D2.^2).*fftnc(b));
  % res = ifftnc((a.D2.^2).*fftnc(b.*a.mask)).*a.mask;
%  res = ifftnc((a.D2.^2).*fftnc(b.*a.mask));
  res = (a.D2).*fftnc(b.*a.mask);
 %res = fftnc(b.*a.mask);
end

