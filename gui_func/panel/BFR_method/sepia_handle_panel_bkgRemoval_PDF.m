%% h = sepia_handle_panel_bkgRemoval_PDF(hParent,h,position)
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
function h = sepia_handle_panel_bkgRemoval_PDF(hParent,h,position)

%% set default values
defaultTol      = 0.1;
defaultMaxIter  = 50;
defaultPadSize  = 40;

%% Tooltips
tooltip.bkgRemoval.PDF.tol      = 'Stopping criteria';
tooltip.bkgRemoval.PDF.maxIter 	= 'Maximum number of iterations allowed';
tooltip.bkgRemoval.PDF.padSize 	= 'No. of zeros to be added after the last array element';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of PDF panel chrildren
h.bkgRemoval.panel.PDF = uipanel(hParent,...
    'Title','Projection onto dipole field (PDF)',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of PDF panel 

    panelParent = h.bkgRemoval.panel.PDF;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % row 1, col 1
    % text|edit field pair: tolerance
    [h.bkgRemoval.PDF.text.tol,h.bkgRemoval.PDF.edit.tol] = sepia_construct_text_edit(...
        panelParent,'Tolerance:',       defaultTol,     [left(1) bottom(1) width height], wratio);
    
    % row 2, col 1
    % text|edit field pair: maximium number of iterations
    [h.bkgRemoval.PDF.text.maxIter,h.bkgRemoval.PDF.edit.maxIter] = sepia_construct_text_edit(...
        panelParent,'Max. iterations:', defaultMaxIter, [left(1) bottom(2) width height], wratio);
   
    % row 3, col 1
    % text|edit field pair: zero padding size
    [h.bkgRemoval.PDF.text.padSize,h.bkgRemoval.PDF.edit.padSize] = sepia_construct_text_edit(...
        panelParent,'Zeropad size:',    defaultPadSize, [left(1) bottom(3) width height], wratio);
    

%% set tooltips
set(h.bkgRemoval.PDF.text.tol,      'Tooltip',tooltip.bkgRemoval.PDF.tol);
set(h.bkgRemoval.PDF.text.maxIter,  'Tooltip',tooltip.bkgRemoval.PDF.maxIter);
set(h.bkgRemoval.PDF.text.padSize,  'Tooltip',tooltip.bkgRemoval.PDF.padSize);

%% Set callbacks
set(h.bkgRemoval.PDF.edit.maxIter,	'Callback', {@EditInputMinMax_Callback,defaultMaxIter,  1,1});
set(h.bkgRemoval.PDF.edit.tol,    	'Callback', {@EditInputMinMax_Callback,defaultTol,      0,0});
set(h.bkgRemoval.PDF.edit.padSize, 	'Callback', {@EditInputMinMax_Callback,defaultPadSize,  1,0});

end