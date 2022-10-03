%% object = getLargestObject(bw,conn)
%
% Input
% --------------
% bw            : binary mask contains objects
% conn          : number of connected component required
%
% Output
% --------------
% object        : mask of the largest object
%
% Description: get the largest connected object from binary mask
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 18 December 2017
% Date last modified:
%
%
function object = getLargestObject(bw,conn)

if nargin<2
    if ndims(bw)>2
        conn = 26;
    else
        conn = 8;
    end
end

bw = bw ~= 0;

CC = bwconncomp(bw,conn);

prevMaxSize = length(CC.PixelIdxList{1,1});
maxComp=1;
for k=1:length(CC.PixelIdxList)
    currSize = length(CC.PixelIdxList{1,k});
    currMaxzSize = max(currSize,prevMaxSize);
    if currMaxzSize~=prevMaxSize
        maxComp = k;
    end
    prevMaxSize = currMaxzSize;
end
object = zeros(size(bw));
object(CC.PixelIdxList{1,maxComp}) = 1;

end