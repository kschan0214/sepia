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
    subwidth(1) = 0.5;
    subwidth(2) = 1-subwidth(1);
    % text|popup pair: select method
    h.bkgRemoval.text.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','text','String','Method:',...
        'units','normalized','position',[left(1) 0.85 width*subwidth(1) height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.BFR.panel.method);
    h.bkgRemoval.popup.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','popup',...
        'String',methodBFRName,...
        'units','normalized','position',[left(1)+width*subwidth(1) 0.85 width*subwidth(2) height]) ;   
    
    % utility function related to background field removal
    h.bkgRemoval.checkbox.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','checkbox',...
        'String','Remove potential B1 residual phase',...
        'units','normalized','position',[left(1) 0.01 0.48 height],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip',tooltip.BFR.panel.polyfit,...
        'value',defaultPolyfit);
    
    % utility function related to erode local field ROI
    h.bkgRemoval.text.imerode = uicontrol('Parent',h.StepsPanel.bkgRemoval ,...
        'Style','text','String','Erode edge voxel(s):',...
        'units','normalized','position',[0.51 0.01 0.19 height],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'Tooltip',tooltip.BFR.panel.erode);
    h.bkgRemoval.edit.imerode = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','edit',...
        'String',num2str(defaultRadius),...
        'units','normalized','position',[0.71 0.01 0.2 height],...
        'backgroundcolor','white');

    
%% create control panel
% define position and size of all method panels
position_child = [0.01 0.15 0.95 0.65];

% construct all method panels
for k = 1:length(function_BFR_method_panel)
    h = feval(function_BFR_method_panel{k},h.StepsPanel.bkgRemoval,h,position_child);
end

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