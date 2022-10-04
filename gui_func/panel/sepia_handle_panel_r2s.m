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
% Date modified: 13 June 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_r2s(hParent,h,position)

% set up method name displayed on GUI
sepia_universal_variables;

% default value

tooltip.R2s.panel.method    = 'Select a R2* mapping method';
%% layout of the panel
nrow        = 5;
rspacing    = 0.01;
ncol        = 2;
cspacing    = 0.01;
[~,~,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

% Set parent of background removal panel
h.StepsPanel.r2sMethod = uipanel(hParent,...
    'Title','R2* mapping',...
    'fontweight', 'bold',...
    'position',[position(1) position(2) 0.95 0.25],...
    'backgroundcolor',get(h.fig,'color'));

%% design of this panel
    
    height = 0.1;
    wratio = 0.5;
    
    % col 1
    % text|popup pair: select method
    [h.r2s.text.r2sMethod,h.r2s.popup.r2sMethod] = sepia_construct_text_popup(...
        h.StepsPanel.r2sMethod,'Method:', methodR2sName, [left(1) 0.85 width height], wratio);
    
%% create control panel
% define position and size of all method panels
position_child = [0.01 0.17 0.95 0.65];

% construct all method panels
for k = 1:length(function_R2s_method_panel)
    h = feval(function_R2s_method_panel{k},h.StepsPanel.r2sMethod,h,position_child);
end

%% set tooltip
set(h.r2s.text.r2sMethod,       'Tooltip',tooltip.R2s.panel.method);
%% set callback function
set(h.r2s.popup.r2sMethod, 'Callback', {@PopupR2sMethod_Callback,h});

end

%% Callback function
% display corresponding background field removal method's panel
function PopupR2sMethod_Callback(source,eventdata,h)

sepia_universal_variables;
                           
% get selected background removal method
method = source.String{source.Value,1} ;

% switch off all panels first
fields = fieldnames(h.r2s.panel);
for kf = 1:length(fields)
    set(h.r2s.panel.(fields{kf}),    'Visible','off');
end

% switch on only target panel
for k = 1:length(methodR2sName)
    if strcmpi(method,methodR2sName{k})
        set(h.r2s.panel.(fields{k}), 'Visible','on');
        break
    end
end

end
