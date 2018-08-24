function [m w] = dwt3D(x, J, af)

% 3-D Discrete Wavelet Transform
%
% USAGE:
%   w = dwt3D(x, stages, af)
% INPUT:
%   x - N1 by N2 by N3 matrix
%       1) Ni all even
%       2) min(Ni) >= 2^(J-1)*length(af)
%   J - number of stages
%   af  - analysis filters
% OUTPUT:
%   w - cell array of wavelet coefficients
% EXAMPLE:
%   [af, sf] = farras;
%   x = rand(128,64,64);
%   w = dwt3D(x,3,af);
%   y = idwt3D(w,3,sf);
%   err = x-y; 
%   max(max(max(abs(err))))
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% for k = 1:J
%     [x w{k}] = afb3D(x, af, af, af);
% end
% w{J+1} = x;
%%
% modified by Shuai 2011.06.29
m = zeros(size(x));
m_size = size(x);
for k = 1:J
    [x w{k}] = afb3D(double(x), af, af, af);
end
w{J+1} = x;
%%
%output for matrix

% hi{1} = LLH;
% hi{2} = LHL;
% hi{3} = LHH;
% hi{4} = HLL;
% hi{5} = HLH;
% hi{6} = HHL;
% hi{7} = HHH;

for k = 1:J
    ha1 = m_size(1)/2;ha2 = m_size(2)/2;ha3 = m_size(3)/2;
    end1 = m_size(1);end2 = m_size(2);end3 = m_size(3);
    %LLH
    m(1:ha1,1:ha2,ha3+1:end3) = w{k}{1};
    %LHL
    m(1:ha1,ha2+1:end2,1:ha3) = w{k}{2};
    %LHH
    m(1:ha1,ha2+1:end2,ha3+1:end3) = w{k}{3};
    
    %HLL
    m(ha1+1:end1,1:ha2,1:ha3) = w{k}{4};
    %HLH
    m(ha1+1:end1,1:ha2,ha3+1:end3) = w{k}{5};
    %HHL
    m(ha1+1:end1,ha2+1:end2,1:ha3) = w{k}{6};
    %HHH
    m(ha1+1:end1,ha2+1:end2,ha3+1:end3) = w{k}{7};
    m_size = m_size/2;
end
m(1:ha1,1:ha2,1:ha3) = w{J+1};
%%
% sbe =1 ;sen = prod(size(x));
% m(sbe:sen)=x;
% for k = J:-1:1  
%     for p =1:7
%         x = w{k}{p};
%         sbe = sen+1;
%         sen = sbe + prod(size(x)) -1 ;
%         m(sbe:sen) = x;
%     end
% end

end
