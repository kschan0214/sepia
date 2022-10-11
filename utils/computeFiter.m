%% function fiter = computeFiter(s,shat,NUM_MAGN)
% s     - measured signal
% shat  - simulated signal
% NUM_GAGN - no. of phase corrupted echoes:
% NUM_MAGN=0 : complex fitting
% NUM_MAGN=length(s) : magnitude fitting
% NUM_MAGN (0,length(s)) : mixed fitting
%
% Description: Compute the fitter for lsqnonlin
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 
% Date last modified:
%
function fiter = computeFiter(s,shat,NUM_MAGN)
if NUM_MAGN == length(s) % Magnitude fitting
    shat1 = abs(shat);
    s1 = abs(s);
    fiter = shat1(:) - s1(:);
elseif NUM_MAGN == 0 % Complex fitting
        fiter2 = shat(:) - s(:);
        fiter2 = [real(fiter2); imag(fiter2)];
%         fiter2 = [real(fiter2), imag(fiter2)];
        fiter = fiter2;
%         fiter = abs(fiter2);
else
    % Compute mixed fitting fit error
    shat1 = abs(shat(1:NUM_MAGN));
    s1 = abs(s(1:NUM_MAGN));
    shat2 = shat(NUM_MAGN+1:end);
    s2 = s(NUM_MAGN+1:end);

    fiter1 = shat1(:) - s1(:);
    fiter2 = shat2(:) - s2(:);
    fiter2 = [real(fiter2); imag(fiter2)];

    fiter = [fiter1;fiter2];
end
fiter = double(fiter);
end