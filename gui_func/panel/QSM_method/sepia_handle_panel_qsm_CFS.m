%% h = sepia_handle_panel_qsm_CFS(hParent,h,position)
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
% Description: This GUI function creates a panel for CFS method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 4 April 2020
%
%
function h = sepia_handle_panel_qsm_CFS(hParent,h,position)

%% set default values
defaultLambda = 0.13;

%% Tooltips
tooltip.qsm.cfs.lambda_text   	= 'Regularisation parameter to control spatial smoothness of QSM';
tooltip.qsm.cfs.lambda_checkbox	= 'Self estimation of lambda based on L-curve approach';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of CFS panel children

h.qsm.panel.cfs = uipanel(hParent,...
    'Title','Closed-form solution with L2-norm regularisation',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');
        
%% Children of CFS panel
    
    % width of each element in a functional column, in normalised unit of
    % the functional column width
    subwidth(1) = width*0.5;
    subwidth(2) = width-subwidth(1);
    
    % row 1
    % text|edit field pair: regulariation parameter
    h.qsm.cfs.text.lambda = uicontrol('Parent',h.qsm.panel.cfs,...
        'Style','text',...
        'String','Lambda:',...
        'units','normalized','position',[left(1) bottom(1) subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.cfs.lambda_text);
    h.qsm.cfs.edit.lambda = uicontrol('Parent',h.qsm.panel.cfs,'Style','edit',...
        'String',num2str(defaultLambda),...
        'units','normalized','position',[left(1)+subwidth(1) bottom(1) subwidth(2) height],...
        'backgroundcolor','white');
    
    % row 2
    % checkbox of self-estimated regularisation
    h.qsm.cfs.checkbox.lambda = uicontrol('Parent',h.qsm.panel.cfs,...
        'Style','checkbox',...
        'String','Self-optimisation by L-curve approach',...
        'units','normalized','position',[left(1) bottom(2) width height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.qsm.cfs.lambda_checkbox);

%% set callbacks
% edit field
set(h.qsm.cfs.edit.lambda,      'Callback',	{@EditInputMinMax_Callback,defaultLambda,0,0});
% checkbox/edit pair
set(h.qsm.cfs.checkbox.lambda,	'Callback', {@CheckboxEditPair_Callback,h.qsm.cfs.edit.lambda,0});

end