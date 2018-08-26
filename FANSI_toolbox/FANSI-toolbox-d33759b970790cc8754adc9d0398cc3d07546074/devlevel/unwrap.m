function [ out ] = unwrap( phase, voxel_size )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


puw = unwrapLaplacian(phase, size(phase) ,voxel_size);

out = phase;

for i = 1:150
    out_old = out;
out = out + 2*pi*round( (puw - out)/(2*pi) );

if sum(abs(out_old(:)-out(:))) < 1
    break;
end

end

end

