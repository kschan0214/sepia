function [slope1, slope2, intrcp]  = fitgrad(g1, g2, in)

Ainv = zeros(length(g1)+1);
Ainv(:,1) = ones(size(Ainv(:,1)));
for j = 2 : size(Ainv,2)
    Ainv(j:end,j) = [1:1:length(Ainv(j:end,j))];
end

g1 = [0;g1];
g2 = [0;g2];
in = [0;in];

slope1 = Ainv*g1;
slope2 = Ainv*g2;
intrcp = Ainv*in;
end