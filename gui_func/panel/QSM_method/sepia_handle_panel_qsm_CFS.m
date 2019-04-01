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
% Date last modified: 
%
%
function h = sepia_handle_panel_qsm_CFS(hParent,h,position)

%% set default values
defaultLambda = 0.13;

%% Parent handle of CFS panel children

h.qsm.panel.cfs = uipanel(hParent,...
    'Title','Closed-form solution with L2-norm regularisation',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');
        

%% Children of CFS panel
    
    % text|edit field pair: regulariation parameter
    h.qsm.cfs.text.lambda = uicontrol('Parent',h.qsm.panel.cfs,...
        'Style','text',...
        'String','Lambda:',...
        'units','normalized','position',[0.01 0.75 0.2 0.2],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Regularisation parameter to control spatial smoothness of QSM');
    h.qsm.cfs.edit.lambda = uicontrol('Parent',h.qsm.panel.cfs,'Style','edit',...
        'String',num2str(defaultLambda),...
        'units','normalized','position',[0.25 0.75 0.2 0.2],...
        'backgroundcolor','white');

    % checkbox of self-estimated regularisation
    h.qsm.cfs.checkbox.lambda = uicontrol('Parent',h.qsm.panel.cfs,...
        'Style','checkbox',...
        'String','Self-optimisation by L-curve approach',...
        'units','normalized','position',[0.01 0.5 1 0.2],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Self estimation of lambda based on L-curve approach');

%% set callbacks
% edit field
set(h.qsm.cfs.edit.lambda,      'Callback',	{@EditInputMinMax_Callback,defaultLambda,0,0});
% checkbox/edit pair
set(h.qsm.cfs.checkbox.lambda,	'Callback', {@CheckboxEditPair_Callback,h.qsm.cfs.edit.lambda,0});

end