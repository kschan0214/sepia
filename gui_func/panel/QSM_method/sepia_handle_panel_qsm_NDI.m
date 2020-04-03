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
% Description: This GUI function creates a panel for iLSQR method
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

    % width of each element in a functional column, in normalised unit of
    % the functional column width
    subwidth(1) = width*0.5;
    subwidth(2) = width-subwidth(1);
    
    % text|edit field pair: tolerance
     h.qsm.NDI.text.tol = uicontrol('Parent',h.qsm.panel.NDI,'Style','text',...
        'String','Tolerance:',...
        'units','normalized','position',[left(1) bottom(1) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.NDI.tol);
    h.qsm.NDI.edit.tol = uicontrol('Parent',h.qsm.panel.NDI,'Style','edit',...
        'String',num2str(defaultTol),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(1) subwidth(2) height],...
        'backgroundcolor','white');

    % text|edit field pair: maximum iterations
    h.qsm.NDI.text.maxIter = uicontrol('Parent',h.qsm.panel.NDI,'Style','text',...
        'String','Max. iterations:',...
        'units','normalized','position',[left(1) bottom(2) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.NDI.maxIter);
    h.qsm.NDI.edit.maxIter = uicontrol('Parent',h.qsm.panel.NDI,'Style','edit',...
        'String',num2str(defaultMaxIter),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(2) subwidth(2) height],...
        'backgroundcolor','white');

    % text|edit field pair: gradient step size
    h.qsm.NDI.text.stepSize = uicontrol('Parent',h.qsm.panel.NDI,'Style','text',...
        'String','Step size:',...
        'units','normalized','position',[left(1) bottom(3) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.NDI.stepSize);
    h.qsm.NDI.edit.stepSize = uicontrol('Parent',h.qsm.panel.NDI,'Style','edit',...
        'String',num2str(defaultStepSize),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(3) subwidth(2) height],...
        'backgroundcolor','white');


%% set callbacks
% edit fields
set(h.qsm.NDI.edit.stepSize,	'Callback',	{@EditInputMinMax_Callback,defaultStepSize,	0,0});
set(h.qsm.NDI.edit.maxIter, 	'Callback',	{@EditInputMinMax_Callback,defaultMaxIter,  1,1});
set(h.qsm.NDI.edit.tol,       	'Callback',	{@EditInputMinMax_Callback,defaultTol,      0,0});

end