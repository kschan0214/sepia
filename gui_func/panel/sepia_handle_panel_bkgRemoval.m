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
% Date last modified: 1 June 2018
%
%
function h = sepia_handle_panel_bkgRemoval(hParent,h,position)
% set up method name displayed on GUI
methodName = {'LBV','PDF','RESHARP','SHARP','VSHARP STI suite','VSHARP','iHARPERELLA'};

% default value
defaultRadius = 0;

% Set parent of background removal panel
h.StepsPanel.bkgRemoval = uipanel(hParent,...
    'Title','Background field removal',...
    'position',[position(1) position(2) 0.95 0.25],...
    'backgroundcolor',get(h.fig,'color'));

%% design of this panel

    % text|popup pair: select method
    h.bkgRemoval.text.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','text','String','Method:',...
        'units','normalized','position',[0.01 0.85 0.3 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Select background field removal method');
    h.bkgRemoval.popup.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','popup',...
        'String',methodName,...
        'units','normalized','position',[0.31 0.85 0.4 0.1]) ;   
    
    % utility function related to background field removal
    h.bkgRemoval.checkbox.bkgRemoval = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','checkbox',...
        'String','Remove B1 residual phase',...
        'units','normalized','position',[0.01 0.01 0.48 0.1],...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Roemer coil combination fails to remove B1+ phase and includes some contribution from B1- of the body coil. These can be removed by fitting a 4-th order 3D polynomial on the computed local field map',...
        'value',1);
    
    % utility function related to erode local field ROI
    h.bkgRemoval.text.imerode = uicontrol('Parent',h.StepsPanel.bkgRemoval ,...
        'Style','text','String','Erode (voxel):',...
        'units','normalized','position',[0.51 0.01 0.19 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'));
    h.bkgRemoval.edit.imerode = uicontrol('Parent',h.StepsPanel.bkgRemoval,...
        'Style','edit',...
        'String',num2str(defaultRadius),...
        'units','normalized','position',[0.71 0.01 0.2 0.1],...
        'backgroundcolor','white',...
        'tooltip','');

    
%% create control panel

% define position and size of all method panels
position_child = [0.01 0.15 0.95 0.65];

    % LBV
    h = sepia_handle_panel_bkgRemoval_LBV(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
        
    % PDF
    h = sepia_handle_panel_bkgRemoval_PDF(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
        
    % SHARP
    h = sepia_handle_panel_bkgRemoval_SHARP(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
    
    % RESHARP    
    h = sepia_handle_panel_bkgRemoval_RESHARP(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
        
    % VSHARPSTI    
    h = sepia_handle_panel_bkgRemoval_VSHARPSTI(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);
    
    % VSHARP
    h = sepia_handle_panel_bkgRemoval_VSHARP(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);

    % iHARPERELLA    
    h = sepia_handle_panel_bkgRemoval_iHARPERELLA(h.StepsPanel.bkgRemoval,...
                                                    h, position_child);

    % in future, add panel of new method here

%% set callback function
set(h.bkgRemoval.popup.bkgRemoval, 'Callback', {@PopupBkgRemoval_Callback,h});
set(h.bkgRemoval.edit.imerode,     'Callback', {@EditInputMinMax_Callback,defaultRadius,1,0});

end

%% Callback function
function PopupBkgRemoval_Callback(source,eventdata,h)
% display corresponding background field removal method's panel

% global h

% get selected background removal method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.bkgRemoval.panel);
for kf = 1:length(fields)
    set(h.bkgRemoval.panel.(fields{kf}),    'Visible','off');
end

% switch on target panel
switch method
    case 'LBV'
        set(h.bkgRemoval.panel.LBV,         'Visible','on');
        
    case 'PDF'
        set(h.bkgRemoval.panel.PDF,         'Visible','on');

    case 'RESHARP'
        set(h.bkgRemoval.panel.RESHARP,     'Visible','on');

    case 'SHARP'
        set(h.bkgRemoval.panel.SHARP,       'Visible','on');

    case 'VSHARP'
        set(h.bkgRemoval.panel.VSHARP,      'Visible','on');

    case 'VSHARP STI suite'
        set(h.bkgRemoval.panel.VSHARPSTI,   'Visible','on');

    case 'iHARPERELLA'
        set(h.bkgRemoval.panel.iHARPERELLA, 'Visible','on');

    % in the future, add new method here
end

end