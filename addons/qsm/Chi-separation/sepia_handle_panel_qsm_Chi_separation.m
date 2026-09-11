%% h = sepia_handle_panel_qsm_Chi_separation(hParent,h,position)
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
% Description: This GUI function creates a panel for Chi separation method
%
% Date created: 02 Dec 2025 by Taechang Kim @ SNU (sakkar2@snu.ac.kr)
% Date modified:
%
%
function h = sepia_handle_panel_qsm_Chi_separation(hParent,h,position)

%% set default values
defaultDr = 137;

menuSolver       = {'Chi-separation-MEDI', 'Chi-separation-iLSQR', 'Chi-sepnet-R2*', 'Chi-sepnet-R2'''};

open_icon = sepia_theme_open_icon();

%% Tooltips
tooltip.qsm.Chi_separation.solver   	 = 'Select a Chi-separation algorithm';
tooltip.qsm.Chi_separation.R2s           = 'Select a R2s map';
tooltip.qsm.Chi_separation.R2            = 'Select a R2 map registered to R2* (GRE)';

%% layout of the panel
nrow        = 4;
rspacing    = 0.03;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent handle of Chi-separation panel children

h.qsm.panel.Chi_separation = uipanel(hParent,...
    'Title','Chi-separation',...
    'position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');

%% Children of Chi-separation panel

    panelParent = h.qsm.panel.Chi_separation;

    % width of each element in a functional column, in normalised unit
    wratio = 0.5;
    
    % col 1, row 1
    [h.qsm.Chi_separation.text.solver,h.qsm.Chi_separation.popup.solver] = sepia_construct_text_popup(...
        panelParent, 'Solver:', menuSolver, [left(1) bottom(1) width height], wratio);

    % col 1, row 2
    [h.qsm.Chi_separation.text.R2s,h.qsm.Chi_separation.edit.R2s,h.qsm.Chi_separation.button.R2s] = sepia_construct_text_edit_button(...
        panelParent, 'R2star:', [], open_icon, [left(1) bottom(2) width height], [0.5 0.45 0.05]);

    % col 1, row 3
    [h.qsm.Chi_separation.text.R2,h.qsm.Chi_separation.edit.R2,h.qsm.Chi_separation.button.R2] = sepia_construct_text_edit_button(...
        panelParent, 'R2:', [], open_icon, [left(1) bottom(3) width height], [0.5 0.45 0.05]);

    % col 1, row 4
    [h.qsm.Chi_separation.text.Dr,h.qsm.Chi_separation.edit.Dr] = sepia_construct_text_edit(...
        panelParent, 'Dr:', defaultDr, [left(1) bottom(4) width height], wratio);
    
%% set tooltips
set(h.qsm.Chi_separation.text.solver,         'Tooltip',tooltip.qsm.Chi_separation.solver);
set(h.qsm.Chi_separation.text.R2s   ,         'Tooltip',tooltip.qsm.Chi_separation.R2s);
set(h.qsm.Chi_separation.text.R2    ,         'Tooltip',tooltip.qsm.Chi_separation.R2);

%% set callbacks
set(h.qsm.Chi_separation.edit.Dr,             'Callback', {@EditInputMinMax_Callback,defaultDr,0,0});
set(h.qsm.Chi_separation.button.R2s,          'Callback', {@ButtonOpen_Callback_,h,'R2star'});
set(h.qsm.Chi_separation.button.R2,           'Callback', {@ButtonOpen_Callback_,h,'R2'});
set(h.qsm.Chi_separation.popup.solver,        'Callback', {@PopupChisepAlgorithm_Callback,h,menuSolver});

end

%% Callback function
function ButtonOpen_Callback_(source,eventdata,h,field)
% get directory and display it on GUI

% global h

switch field
    case 'R2star'
        % only read NIfTI file for R2star
        [fileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select a NIfTI file for R2star');

        if pathDir ~= 0
            % set input edit field for display
            set(h.qsm.Chi_separation.edit.R2s,    'String',fullfile(pathDir,fileName));
        end

    case 'R2'
        % only read NIfTI file for R2star
        [fileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select a NIfTI file for R2');

        if pathDir ~= 0
            % set input edit field for display
            set(h.qsm.Chi_separation.edit.R2,    'String',fullfile(pathDir,fileName));
        end        

end

end

% Change default Dr value depending on chi-separation algorithm
function PopupChisepAlgorithm_Callback(source,eventdata,h,menuSolver)

% needs 'methodQSMName' here
sepia_universal_variables;

switch source.String{source.Value}
    case menuSolver{1}
        set(h.qsm.Chi_separation.edit.Dr,'String',num2str(137));
        set(h.qsm.Chi_separation.edit.Dr, 'Enable','on');
        set(h.qsm.Chi_separation.edit.R2s,'Enable','on');
        set(h.qsm.Chi_separation.edit.R2, 'Enable','on');
    case menuSolver{2}
        set(h.qsm.Chi_separation.edit.Dr,'String',num2str(137));
        set(h.qsm.Chi_separation.edit.Dr, 'Enable','on');
        set(h.qsm.Chi_separation.edit.R2s,'Enable','on');
        set(h.qsm.Chi_separation.edit.R2, 'Enable','on');
    case menuSolver{3}
        set(h.qsm.Chi_separation.edit.Dr, 'String',num2str(114));
        set(h.qsm.Chi_separation.edit.Dr, 'Enable','off');
        set(h.qsm.Chi_separation.edit.R2s,'Enable','on');
        set(h.qsm.Chi_separation.edit.R2, 'Enable','off');
    case menuSolver{4}
        set(h.qsm.Chi_separation.edit.Dr, 'String',num2str(114));
        set(h.qsm.Chi_separation.edit.Dr, 'Enable','off');
        set(h.qsm.Chi_separation.edit.R2s,'Enable','on');
        set(h.qsm.Chi_separation.edit.R2, 'Enable','on');
end

end