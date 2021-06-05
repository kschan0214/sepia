%% [height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing)
%
% Input
% --------------
% nrow          : number of rows in the panel
% rspacing      : spacing between rows, in normalised unit
% ncol          : number of functional columns in the panel
% cspacing      : spacing between column, in normalised unit
%
% Output
% --------------
% height        : height of a row, in normalised unit
% bottom        : starting positions of rows
% width         : width of a functional column, in normalised unit
% left          : starting positions of function columns
%
% Description: compute the normalised position of GUI elements
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 3 April 2020
% Date modified:
%
%
function [height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing)

if nargin < 4
    cspacing = 0;
end
if nargin < 3
    ncol = 1;
end
if nargin < 2
    rspacing = 0;
end

height  = (1-(nrow+1)*rspacing)/nrow;    % height of a row in relative unit
bottom  = fliplr((height+rspacing:height+rspacing:(height+rspacing)*nrow) - height);     % start positions of rows
width   = (1-(ncol+1)*cspacing)/ncol;       % 
left    = (cspacing:width+cspacing:(width+cspacing)*ncol);  % start positions of columns


end