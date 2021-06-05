%% h = sepia_handle_panel_bkgRemoval_LBV(hParent,h,position)
%
% Input
% --------------
% hParent       : parent handle of this panel
% h             : global structure contains all handles
% position      : position of this panel
%
% Output
% --------------
% h             : global structure contains all new and other handles
%
% Description: This GUI function creates a panel for LBV method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 3 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_bkgRemoval_LBV(hParent,h,position)

%% set default values
defaultTol      = 0.0001;
defaultDepth    = 5;
defaultPeel     = 2;

%% Tooltips
tooltip.bkgRemoval.LBV.tol      = 'Iteration stopping criteria on the coarsest grid';
tooltip.bkgRemoval.LBV.depth  	= 'Multigrid level: no. of length scales. The largest length scale is 2^depth * voxel size';
tooltip.bkgRemoval.LBV.peel     = 'No. of boundary layers to be peeled off (i.e. disgarded)';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of LBV panel children

h.bkgRemoval.panel.LBV = uipanel(hParent,...
    'Title','Laplacian boundary value (LBV)',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'));

%% Children of LBV panel

    panelParent = h.bkgRemoval.panel.LBV;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % functional col: 'text|edit' pair: tolerance
    [h.bkgRemoval.LBV.text.tol,h.bkgRemoval.LBV.edit.tol] = sepia_construct_text_edit(...
        panelParent,'Tolerance:',   defaultTol,     [left(1) bottom(1) width height], wratio);

    % row 2, col 1
    % functional col: 'text|edit' pair: depth
    [h.bkgRemoval.LBV.text.depth,h.bkgRemoval.LBV.edit.depth] = sepia_construct_text_edit(...
        panelParent,'Depth:',       defaultDepth,   [left(1) bottom(2) width height], wratio);
    
    % row 3, col 1
    % functional col: 'text|edit' pair: peel
    [h.bkgRemoval.LBV.text.peel,h.bkgRemoval.LBV.edit.peel] = sepia_construct_text_edit(...
        panelParent,'Peel:',        defaultPeel,    [left(1) bottom(3) width height], wratio);
    
    
%% set tooltips
set(h.bkgRemoval.LBV.text.tol,  'Tooltip',tooltip.bkgRemoval.LBV.tol);
set(h.bkgRemoval.LBV.text.depth,'Tooltip',tooltip.bkgRemoval.LBV.depth);
set(h.bkgRemoval.LBV.text.peel, 'Tooltip',tooltip.bkgRemoval.LBV.peel);
    
%% Set callbacks
set(h.bkgRemoval.LBV.edit.depth,    'Callback', {@EditInputMinMax_Callback,defaultDepth ,1,-1});
set(h.bkgRemoval.LBV.edit.peel,     'Callback', {@EditInputMinMax_Callback,defaultPeel  ,1,0});
set(h.bkgRemoval.LBV.edit.tol,      'Callback', {@EditInputMinMax_Callback,defaultTol   ,0,0});

end