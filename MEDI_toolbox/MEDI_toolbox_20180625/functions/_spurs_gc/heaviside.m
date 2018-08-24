function Y = heaviside(X)
%HEAVISIDE    Step function.
%    HEAVISIDE(X) is 0 for X < 0, 1 for X > 0, and .5 for X == 0.
%    HEAVISIDE(X) is not a function in the strict sense.
%    See also DIRAC.

%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2009/03/09 20:41:30 $

Y = zeros(size(X));
Y(X > 0) = 1;

Y(X == 0) = .5;



% Y = zeros(size(X));
% Y(X > 0) = 1;
% eng = symengine;
% if strcmp(eng.kind,'maple')
%    Y(X == 0) = nan;
% else
%    Y(X == 0) = .5;
% end