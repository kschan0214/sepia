function [ rmse ] = rmse( in, true, use_abs )
%RMSE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    rmse = 100 * norm(in(:) - true(:)) / norm(true(:));
end

if (nargin == 3) && (use_abs == 1)
    rmse = 100 * norm(abs(in(:)) - abs(true(:))) / norm(abs(true(:)));
end

end

