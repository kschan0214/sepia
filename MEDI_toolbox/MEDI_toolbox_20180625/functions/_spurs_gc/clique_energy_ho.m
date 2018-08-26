
function e=clique_energy_ho(d,p,th,quant)
%clique_energy_ho  Computes clique energy: e=th^(p-2)*d.^2.*mask + d.^p.*(1-mask)
%        e=clique_energy_ho(d,p,th,quant)
%
%  Input arguments --------------------
%  d             -> clique difference
%  p             -> power law exponent
%  th            -> it defines a region over which the potential grows quadratically
%  quant         -> it defines whether or not the potential is quantized

switch quant
    case 'no'        % non quantized potential
        d=abs(d);
    case 'yes'       % quantized potential (2pi Quantization of phase difference)
        d=abs(round(d/2/pi)*2*pi);
end;
%

if th~=0
    mask = (d<=1*th);
    
%     e = th^(p-2)*d.^2.*mask + d.^p.*(1-mask);     % l2

%     e = th^(p-2)*d.^2.*mask + d.^p.*(1-mask);
%     e = th^(p-2)*d.^2.*mask + (th*ones(size(d))).^p.*(1-mask);  %  l2 tru
%    e = 0*ones(size(d)).*mask + (th*ones(size(d))).^p.*(1-mask); % this line is 101 potential
     e = d.^p;
 
%     e = d.^p.*mask + (th*ones(size(d))).^p.*(1-mask);  % l1_truncted

else
    e = d.^p;
end

return

