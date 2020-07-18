%% h = sepia_handle_panel_bkgRemoval_RESHARP(hParent,h,position)
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
% Description: This GUI function creates a panel for RESHARP method
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 1 June 2018
% Date modified: 4 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_bkgRemoval_RESHARP(hParent,h,position)

%% set default values
defaultRadius = 4;
defaultLambda = 0.01;

%% Tooltips
tooltip.bkgRemoval.RESHARP.radius  	= 'Radius of spherical mean value kernel';
tooltip.bkgRemoval.RESHARP.lambda  	= 'Regularizaiton parameter used in Tikhonov';

%% layout of the panel
nrow        = 4;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of RESHARP panel children

h.bkgRemoval.panel.RESHARP = uipanel(hParent,...
    'Title','Regularisation SHARP',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of RESHARP panel
    
    panelParent = h.bkgRemoval.panel.RESHARP;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;

    % row 1, col 1
    % text|edit field pair: radius
    [h.bkgRemoval.RESHARP.text.radius,h.bkgRemoval.RESHARP.edit.radius] = sepia_construct_text_edit(...
        panelParent,'SMV radius (voxel):',          defaultRadius,	[left(1) bottom(1) width height], wratio);

    % row 2, col 1
    % text|edit field pair: regularisation parameter
    [h.bkgRemoval.RESHARP.text.lambda,h.bkgRemoval.RESHARP.edit.lambda] = sepia_construct_text_edit(...
        panelParent,'Regularisation parameter:',	defaultLambda,	[left(1) bottom(2) width height], wratio);
    

%% set tooltips
set(h.bkgRemoval.RESHARP.text.radius,	'Tooltip',tooltip.bkgRemoval.RESHARP.radius);
set(h.bkgRemoval.RESHARP.text.lambda,   'Tooltip',tooltip.bkgRemoval.RESHARP.lambda);

%% set callbacks
set(h.bkgRemoval.RESHARP.edit.lambda,	'Callback', {@EditInputMinMax_Callback,defaultLambda,0,0});
set(h.bkgRemoval.RESHARP.edit.radius, 	'Callback', {@EditInputMinMax_Callback,defaultRadius,1,0});

end