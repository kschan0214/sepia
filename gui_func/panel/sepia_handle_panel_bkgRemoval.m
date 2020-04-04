%% h = sepia_handle_panel_bkgRemoval(hParent,h,position)
%
% Input
% --------------
% hParent       : parent handle of this panel
% hFig          : handle of the GUI
% h             : global structure contains all handles
% position      : position of this panel
%
% Output
% --------------
% h             : global structure contains all new and other handles
%
% Description: This GUI function creates a panel for background field
% removal method control
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 April 2018
% Date modified: 1 June 2018
% Date modified: 3 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_bkgRemoval(hParent,h,position)

% set up method name displayed on GUI
sepia_universal_variables;

% default value
defaultRadius   = 0;
defaultPolyfit  = true;

tooltip.BFR.panel.method    = 'Select a background field removal method';
tooltip.BFR.panel.polyfit   = 'Remove a 4th order 3D polynomial from the local field map which might be contributed from residual B1 phase.';
tooltip.BFR.panel.erode     = 'Remove edge voxels from the local field. Might improve the final QSM result';

%% layout of the panel
nrow        = 5;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[~,~,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%%
% Set parent of background removal panel
h.StepsPanel.bkgRemoval = uipanel(hParent,...
    'Title','Background field removal',...
    'fontweight', 'bold',...
    'position',[position(1) position(2) 0.95 0.25],...
    'backgroundcolor',get(h.fig,'color'));

%% design of this panel
    
    height = 0.1;
    wratio = 0.5;
    
    % col 1
    % text|popup pair: select method
    [h.bkgRemoval.text.bkgRemoval,h.bkgRemoval.popup.bkgRemoval] = sepia_construct_text_popup(...
        h.StepsPanel.bkgRemoval,'Method:', methodBFRName, [left(1) 0.85 width height], wratio);

    % utility function related to background field removal
    h.bkgRemoval.checkbox.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','checkbox',...
        'String','Remove potential B1 residual phase',...
        'units','normalized','position',[left(1) 0.01 width height],...
        'backgroundcolor',get(h.fig,'color'),...
        'value',defaultPolyfit);
    
    % col 2
    % text|field pair: utility function related to erode local field ROI
    [h.bkgRemoval.text.imerode,h.bkgRemoval.edit.imerode] = sepia_construct_text_edit(...
        h.StepsPanel.bkgRemoval,'Erode edge voxel(s):', defaultRadius, [left(2) 0.01 width height], wratio);

    
%% create control panel
% define position and size of all method panels
position_child = [0.01 0.15 0.95 0.65];

% construct all method panels
for k = 1:length(function_BFR_method_panel)
    h = feval(function_BFR_method_panel{k},h.StepsPanel.bkgRemoval,h,position_child);
end

%% set tooltip
set(h.bkgRemoval.text.bkgRemoval,       'Tooltip',tooltip.BFR.panel.method);
set(h.bkgRemoval.checkbox.bkgRemoval,	'Tooltip',tooltip.BFR.panel.polyfit);
set(h.bkgRemoval.text.imerode,          'Tooltip',tooltip.BFR.panel.erode);

%% set callback function
set(h.bkgRemoval.popup.bkgRemoval, 'Callback', {@PopupBkgRemoval_Callback,h});
set(h.bkgRemoval.edit.imerode,     'Callback', {@EditInputMinMax_Callback,defaultRadius,1,0});

end

%% Callback function
% display corresponding background field removal method's panel
function PopupBkgRemoval_Callback(source,eventdata,h)

sepia_universal_variables;
                           
% get selected background removal method
method = source.String{source.Value,1} ;

% switch off all panels first
fields = fieldnames(h.bkgRemoval.panel);
for kf = 1:length(fields)
    set(h.bkgRemoval.panel.(fields{kf}),    'Visible','off');
end

% switch on only target panel
for k = 1:length(methodBFRName)
    if strcmpi(method,methodBFRName{k})
        set(h.bkgRemoval.panel.(fields{k}), 'Visible','on');
        break
    end
end

end