%% h = sepia_handle_panel_qsm_NDI(hParent,h,position)
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
% Description: This GUI function creates a panel for NDI method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 5 June 2019
% Date modified: 4 April 2020
%
%
function h = sepia_handle_panel_qsm_NDI(hParent,h,position)

%% set default values
defaultTol      = 1;
defaultMaxIter  = 200;
defaultStepSize = 1;

%% Tooltips
tooltip.qsm.NDI.tol   	 = 'Relative tolerance to stop NDI iteration.';
tooltip.qsm.NDI.maxIter	 = 'Maximum iterations allowed';
tooltip.qsm.NDI.stepSize = 'Gradient descent step size';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of CFS panel children

h.qsm.panel.NDI = uipanel(hParent,...
    'Title','Non-linear Dipole Inversion (NDI)',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of CFS panel

    panelParent = h.qsm.panel.NDI;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % col 1, row 1
    % text|edit field pair: tolerance
    [h.qsm.NDI.text.tol,h.qsm.NDI.edit.tol] = sepia_construct_text_edit(...
        panelParent,'Tolerance:',       defaultTol,      [left(1) bottom(1) width height], wratio);
    
    % col 1, row 2
    % text|edit field pair: maximum iterations
    [h.qsm.NDI.text.maxIter,h.qsm.NDI.edit.maxIter] = sepia_construct_text_edit(...
        panelParent,'Max. iterations:', defaultMaxIter,  [left(1) bottom(2) width height], wratio);
    
    % col 1, row 3
    % text|edit field pair: gradient step size
    [h.qsm.NDI.text.stepSize,h.qsm.NDI.edit.stepSize] = sepia_construct_text_edit(...
        panelParent,'Step size:',       defaultStepSize, [left(1) bottom(3) width height], wratio);
    
    
%% set tooltips
set(h.qsm.NDI.text.tol,         'Tooltip',tooltip.qsm.NDI.tol);
set(h.qsm.NDI.text.maxIter,     'Tooltip',tooltip.qsm.NDI.maxIter);
set(h.qsm.NDI.text.stepSize,	'Tooltip',tooltip.qsm.NDI.stepSize);

%% set callbacks
% edit fields
set(h.qsm.NDI.edit.stepSize,	'Callback',	{@EditInputMinMax_Callback,defaultStepSize,	0,0});
set(h.qsm.NDI.edit.maxIter, 	'Callback',	{@EditInputMinMax_Callback,defaultMaxIter,  1,1});
set(h.qsm.NDI.edit.tol,       	'Callback',	{@EditInputMinMax_Callback,defaultTol,      0,0});

end