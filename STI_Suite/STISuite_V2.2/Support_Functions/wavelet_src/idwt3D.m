function y = idwt3D(J, sf , m ,w)

% Inverse 3-D Discrete Wavelet Transform
%
% USAGE:
%   y = idwt3D(w, J, sf)
% INPUT:
%   w - wavelet coefficient
%   J  - number of stages
%   sf - synthesis filters
% OUTPUT:
%   y - output array
% See: dwt3D
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/
% modified by Shuai 2011.06.29
% y = w{J+1};
% for k = J:-1:1
%    y = sfb3D(y, w{k}, sf, sf, sf);
% end
%%
m_size = size(m);
if (nargin<4)
    for k = 1:J
        ha1 = m_size(1)/2;ha2 = m_size(2)/2;ha3 = m_size(3)/2;
        end1 = m_size(1);end2 = m_size(2);end3 = m_size(3);
        %LLH
        w{k}{1}= m(1:ha1,1:ha2,ha3+1:end3) ;
        %LHL
        w{k}{2}= m(1:ha1,ha2+1:end2,1:ha3) ;
        %LHH
        w{k}{3}= m(1:ha1,ha2+1:end2,ha3+1:end3) ;
        
        %HLL
        w{k}{4}= m(ha1+1:end1,1:ha2,1:ha3) ;
        %HLH
        w{k}{5}= m(ha1+1:end1,1:ha2,ha3+1:end3) ;
        %HHL
        w{k}{6}= m(ha1+1:end1,ha2+1:end2,1:ha3) ;
        %HHH
        w{k}{7}= m(ha1+1:end1,ha2+1:end2,ha3+1:end3) ;
        m_size = m_size/2;
    end
    w{J+1} = m(1:ha1,1:ha2,1:ha3);
end
y = w{J+1};
for k = J:-1:1
    y = sfb3D(y, w{k}, sf, sf, sf);
end
end

