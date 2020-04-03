%% h = sepia_handle_panel_qsm_iLSQR(hParent,h,position)
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
% Date created: 1 June 2018
% Date modified: 4 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_qsm_iLSQR(hParent,h,position)

%% set default values
defaultTol      = 0.001;
defaultLambda   = 0.13;
defaultMaxIter  = 100;

%% Tooltips
tooltip.qsm.iLSQR.tol               = 'Relative tolerance to stop LSQR iteration';
tooltip.qsm.iLSQR.maxIter           = 'Maximum iterations allowed';
tooltip.qsm.iLSQR.lambda_text       = 'Regularisation parameter to control spatial smoothness of QSM';
tooltip.qsm.iLSQR.lambda_checkbox	= 'Self estimation of lambda based on L-curve approach';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of CFS panel children

h.qsm.panel.iLSQR = uipanel(hParent,...
    'Title','Iterative LSQR',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of CFS panel

    % width of each element in a functional column, in normalised unit of
    % the functional column width
    subwidth(1) = width*0.5;
    subwidth(2) = width-subwidth(1);
    
    % functional colume 1
    % text|edit field pair: tolerance
     h.qsm.iLSQR.text.tol = uicontrol('Parent',h.qsm.panel.iLSQR,'Style','text',...
        'String','Tolerance:',...
        'units','normalized','position',[left(1) bottom(1) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.iLSQR.tol);
    h.qsm.iLSQR.edit.tol = uicontrol('Parent',h.qsm.panel.iLSQR,'Style','edit',...
        'String',num2str(defaultTol),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(1) subwidth(2) height],...
        'backgroundcolor','white');

    % text|edit field pair: maximum iterations
    h.qsm.iLSQR.text.maxIter = uicontrol('Parent',h.qsm.panel.iLSQR,'Style','text',...
        'String','Max. iterations:',...
        'units','normalized','position',[left(1) bottom(2) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.iLSQR.maxIter);
    h.qsm.iLSQR.edit.maxIter = uicontrol('Parent',h.qsm.panel.iLSQR,'Style','edit',...
        'String',num2str(defaultMaxIter),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(2) subwidth(2) height],...
        'backgroundcolor','white');

    % text|edit field pair: regularisation parameter
    h.qsm.iLSQR.text.lambda = uicontrol('Parent',h.qsm.panel.iLSQR,'Style','text',...
        'String','Lambda:',...
        'units','normalized','position',[left(1) bottom(3) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.iLSQR.lambda_text);
    h.qsm.iLSQR.edit.lambda = uicontrol('Parent',h.qsm.panel.iLSQR,'Style','edit',...
        'String',num2str(defaultLambda),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(3) subwidth(2) height],...
        'backgroundcolor','white');

    % checkbox of self-estimated regularisation
    h.qsm.iLSQR.checkbox.lambda = uicontrol('Parent',h.qsm.panel.iLSQR,'Style','checkbox',...
        'String','Self-optimisation by L-curve approach',...
        'units','normalized','position',[left(1) bottom(4) width height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.iLSQR.lambda_checkbox);

%% set callbacks
% edit fields
set(h.qsm.iLSQR.edit.lambda,        'Callback',	{@EditInputMinMax_Callback,defaultLambda,   0,0});
set(h.qsm.iLSQR.edit.maxIter,       'Callback',	{@EditInputMinMax_Callback,defaultMaxIter,  1,1});
set(h.qsm.iLSQR.edit.tol,           'Callback',	{@EditInputMinMax_Callback,defaultTol,      0,0});

% checkbox/edit pair
set(h.qsm.iLSQR.checkbox.lambda,	'Callback', {@CheckboxEditPair_Callback,h.qsm.iLSQR.edit.lambda,0});

end