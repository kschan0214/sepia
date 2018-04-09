function x = L2_Linf_shrink(y,t)

% This function minimizes
%     0.5*||b*x-y||_2^2 + t*||x||_inf
% where b is a scalar.  Note that it suffices
% to consider the minimization
%     0.5*||x-y||_2^2 + t/b*||x||_inf
% and so we will assume that the value of b has
% been absorbed into t (= tau).
% The minimization proceeds by initializing
% x with y.  Let z be y re-ordered so that
% the abs(z) is in descending order.  Then
% first solve
%     min_{b>=abs(z2)} 0.5*(b-abs(z1))^2 + t*b
% if b* = abs(z2), then repeat with first and
% second largest z values;
%     min_{b>=abs(z3)} 0.5*(b-abs(z1))^2+0.5*(b-abs(z2))^2 + t*b
% which by expanding the square is equivalent to
%     min_{b>=abs(z3)} 0.5*(b-mean(abs(z1),abs(z2)))^2 + t*b
% and repeat this process if b*=abs(z3), etc.
% This reduces problem to finding a cut-off index, where
% all coordinates are shrunk up to and including that of
% the cut-off index.  The cut-off index is the smallest
% integer k such that
%    1/k sum(abs(z(1)),...,abs(z(k))) - t/k <= abs(z(k+1))
%

x = y;
[dummy,o] = sort(abs(y),'descend');
z = y(o);
mz = abs(z);

% find cut-off index
cs = cumsum(abs(z(1:length(z)-1)))./(1:length(z)-1)'-t./(1:length(z)-1)';
d = (cs>abs(z(2:length(z))));
if sum(d) == 0
   cut_index = length(y);
else
   cut_index = min(find(d==1));
end

% shrink coordinates 1 to cut_index
zbar = mean(abs(z(1:cut_index)));
if cut_index < length(y)
   x(o(1:cut_index)) = sign(z(1:cut_index))*max(zbar-t/cut_index,abs(z(cut_index+1)));
else
   x(o(1:cut_index)) = sign(z(1:cut_index))*max(zbar-t/cut_index,0);
end
