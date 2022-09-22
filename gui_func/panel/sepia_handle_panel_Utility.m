%% h = sepia_handle_panel_Utility(hParent,h,position)
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
function h = sepia_handle_panel_Utility(hParent,h,position)

% set up method name displayed on GUI
utilityName = {'Get header info','Get lateral ventricle mask','Manage Dependency','Convert GE real/imaginary images to phase image (experimental)'};

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
    'position',[position(1) position(2) 0.98 0.60]);
    
%     % Temporo-spatial unwrapping methods
%     h.utility.text.utilityMethod = uicontrol('Parent',h.StepsPanel.Utility,...
%         'Style','text','String','Utility:',...
%         'units','normalized','position',[0.01 0.8 0.3 0.2],...
%         'HorizontalAlignment','left',...
%         'backgroundcolor',get(h.fig,'color'),...
%         'tooltip','Select utility method');
    h.Utility.popup.utilityMethod = uicontrol('Parent',h.StepsPanel.Utility,...
        'Style','popup',...
        'String',utilityName,...
        'units','normalized','position',[0.31 0.85 0.5 0.1]); 
    
    
    % define position and size of all method panels
position_child = [0.01 0.05 0.98 0.75];

    % get header
    h = sepia_handle_panel_utility_get_header(h.StepsPanel.Utility,...
                                                    h,position_child);
                                                
    h = sepia_handle_panel_utility_mask_ventricle(h.StepsPanel.Utility,...
                                                    h,position_child);

    h = sepia_handle_panel_utility_manage_dependency(h.StepsPanel.Utility,...
                                                    h,position_child);
                                                
    h = sepia_handle_panel_utility_convert_realImaginary2phase(h.StepsPanel.Utility,...
                                                    h,position_child);


% %% set callback functions
set(h.Utility.popup.utilityMethod, 'Callback', {@popupUtility_Callback,h});
        
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
        set(h.Utility.panel.getHeader,          'Visible','on');
        
    case 'Get lateral ventricle mask'
        set(h.Utility.panel.csfMask,            'Visible','on');
        
    case 'Manage Dependency'
        set(h.Utility.panel.magageDependency,	'Visible','on');
        
    case 'Convert GE real/imaginary images to phase image (experimental)'
        set(h.Utility.panel.realimag2phase,     'Visible','on');

    % in the future, add new method here
end

end