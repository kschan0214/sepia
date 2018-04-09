% Aswitch.m
%
% Add switch for spgl1
%

function y = Aswitch(x, mode, A, At)

if (mode == 1)
  y = A(x);
elseif (mode == 2)
  y = At(x);
else
  error('There is a problem.');
end
