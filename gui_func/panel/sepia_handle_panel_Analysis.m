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
% Date created: 9 October 2022
% Date modified:
%
%
function h = sepia_handle_panel_Analysis(hParent,h,position)
global analysisName
% set up method name displayed on GUI
analysisName = {'Segmentation - CIT168 Reinf. learn atlas',...
                'Segmentation - MuSus100 atlas',...
                'Segmentation - AHEAD'};

% Parent handle of phase unwrapping panel
h.StepsPanel.Analysis = uipanel(hParent,'Title','Analysis',...
    'position',[position(1) position(2) 0.98 0.98]);
    
    h.Analysis.popup.analysisMethod = uicontrol('Parent',h.StepsPanel.Analysis,...
        'Style','popup',...
        'String',analysisName,...
        'units','normalized','position',[0.31 0.94 0.5 0.05]); 
    
    
    % define position and size of all method panels
position_child = [0.01 0.01 0.98 0.9];

    % CIT168 rei. lea. atlas
    h = sepia_handle_panel_analysis_segmentation_CIT168RL(h.StepsPanel.Analysis,...
                                                    h,position_child);
                                                
    % MuSus100
    h = sepia_handle_panel_analysis_segmentation_MuSus100(h.StepsPanel.Analysis,...
                                                    h,position_child);
                                                
    % AHEAD
    h = sepia_handle_panel_analysis_segmentation_AHEAD(h.StepsPanel.Analysis,...
                                                    h,position_child);

%% set callback functions
set(h.Analysis.popup.analysisMethod, 'Callback', {@popupAnalysis_Callback,h});
        
end

function popupAnalysis_Callback(source,eventdata,h)
global analysisName

% get selected background removal method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.Analysis.panel);
for kf = 1:length(fields)
    set(h.Analysis.panel.(fields{kf}),    'Visible','off');
end

% switch on target panel
switch method
    case analysisName{1}
        set(h.Analysis.panel.Segmentation_CIT168RL,          'Visible','on');
        
    case analysisName{2}
        set(h.Analysis.panel.Segmentation_MuSus100,          'Visible','on');
        
    case analysisName{3}
        set(h.Analysis.panel.Segmentation_AHEAD,          'Visible','on');
    % in the future, add new method here
end

end