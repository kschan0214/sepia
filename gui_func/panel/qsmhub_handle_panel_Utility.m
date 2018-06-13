%% h = qsmhub_handle_panel_Utility(hParent,h,position)
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
% Description: This GUI function creates a panel for phase unwrapping method control
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 12 June 2018
% Date modified:
%
%
function h = qsmhub_handle_panel_Utility(hParent,h,position)

% set up method name displayed on GUI
utilityName = {'Get header info','Get lateral ventricle mask'};

% % set default value
% defaultThreshold = 0.5;
% 
% % define maximum level of options and spacing between options
% nlevel = 5;
% spacing = 0.02;
% height = (1-(nlevel+1)*spacing)/nlevel;
% button = (height+spacing:height+spacing:(height+spacing)*nlevel) - height;

% Parent handle of phase unwrapping panel
h.StepsPanel.Utility = uipanel(hParent,'Title','Utility',...
    'position',[position(1) position(2) 0.98 0.70]);
    
%     % Temporo-spatial unwrapping methods
%     h.utility.text.utilityMethod = uicontrol('Parent',h.StepsPanel.Utility,...
%         'Style','text','String','Utility:',...
%         'units','normalized','position',[0.01 0.8 0.3 0.2],...
%         'HorizontalAlignment','left',...
%         'backgroundcolor',get(h.fig,'color'),...
%         'tooltip','Select utility method');
    h.utility.popup.utilityMethod = uicontrol('Parent',h.StepsPanel.Utility,...
        'Style','popup',...
        'String',utilityName,...
        'units','normalized','position',[0.31 0.85 0.5 0.1]); 
    
    
    % define position and size of all method panels
position_child = [0.01 0.24 0.98 0.6];

    % get header
    h = qsmhub_handle_panel_utility_get_header(h.StepsPanel.Utility,...
                                                    h,position_child);
                                                
    h = qsmhub_handle_panel_utility_mask_ventricle(h.StepsPanel.Utility,...
                                                    h,position_child);


% %% set callback functions
set(h.utility.popup.utilityMethod, 'Callback', {@popupUtility_Callback,h});
        
end

function popupUtility_Callback(source,eventdata,h)

% get selected background removal method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.Utility.panel);
for kf = 1:length(fields)
    set(h.Utility.panel.(fields{kf}),    'Visible','off');
end

% switch on target panel
switch method
    case 'Get header info'
        set(h.Utility.panel.getHeader,         'Visible','on');
        
    case 'Get lateral ventricle mask'
        set(h.Utility.panel.csfMask,         'Visible','on');

    % in the future, add new method here
end

end