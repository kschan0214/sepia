%% h = sepia_handle_panel_dataIO(hParent,h,position)
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
% Description: This GUI function creates a panel for I/O control
%
% Kwok-Shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 16 April 2018
% Date modified: 30 May 2018
% Date modified: 5 April 2020 (v0.8.0)
%
%
function h = sepia_handle_panel_dataIO(hParent,h,position)

% default value
defaultFractionalThres  = 0.5;
defaultGradientThres    = 0; 

open_icon = imread('folder@0,3x.jpg');
open_icon = imresize(open_icon,[1 1]*16);

%% Tooltips
tooltip.dataIO.SEPIA.input   	= 'Directory contains all essential files (*ph*.nii*, *mag*.nii* and *header*.mat)';
tooltip.dataIO.SEPIA.output     = 'Output directory with filename prefix. Default output prefix is ''/input_dir/output/sepia''';
tooltip.dataIO.SEPIA.maskdir	= ['Specify a NIfTI mask file (.nii/.nii.gz). ',...
                                   'If the input directory contains a NIfTI file with string *mask* in the filename then it will be read automatically'];
tooltip.dataIO.SEPIA.inputData1	= 'Specify a NIfTI phase image (.nii/.nii.gz)';
tooltip.dataIO.SEPIA.inputData2	= 'Specify a NIfTI magnitude image (.nii/.nii.gz)';
tooltip.dataIO.SEPIA.inputData3	= 'Specify a NIfTI weight map for QSM algorithm (.nii/.nii.gz)';
tooltip.dataIO.SEPIA.inputHeader= 'Specify a SEPIA header';
tooltip.dataIO.SEPIA.invertPhase= 'Reverse the direction of frequency shift. This will invert the contrast of the output.';
tooltip.dataIO.SEPIA.BET        = 'Use FSL bet for brain extraction';
tooltip.dataIO.SEPIA.fractThres = 'Fractional intensity threshold (0->1); default=0.5; smaller values give larger brain outline estimates';
tooltip.dataIO.SEPIA.gradThres  = 'Vertical gradient in fractional intensity threshold (-1->1); default=0; positive values give larger brain outline at bottom, smaller at top';

%% layout of the panel
nrow        = 5;
rspacing    = 0.02;
ncol        = 2;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%% Parent of dataIO panel
h.StepsPanel.dataIO = uipanel(hParent,'Title','I/O',...
    'Position',[position(1) position(2) 0.95 0.2],...
    'fontweight', 'bold',...
    'backgroundcolor',get(h.fig,'color'));

    panelParent = h.StepsPanel.dataIO;
    
    wratio = [0.3,0.65,0.05];
    
    % row 1, col 1
    % input directory
    pos = [left(1) bottom(1) width height];
    [h.dataIO.text.input,h.dataIO.edit.input,h.dataIO.button.input] = sepia_construct_text_edit_button(panelParent,...
        'Input directory:',[],open_icon,pos,wratio);

    % row 1, col 2
    % output basename
    pos = [left(2) bottom(1) width height];
    [h.dataIO.text.output,h.dataIO.edit.output,h.dataIO.button.output] = sepia_construct_text_edit_button(panelParent,...
        'Output prefix:',[],open_icon,pos,wratio);

    % row 2, col 2
    % brain mask
    pos = [left(2) bottom(2) width height];
    [h.dataIO.text.maskdir,h.dataIO.edit.maskdir,h.dataIO.button.maskdir] = sepia_construct_text_edit_button(panelParent,...
        'Brain mask:',[],open_icon,pos,wratio);

    % row 2, col 1
    % input data1
    pos = [left(1) bottom(2) width height];
    [h.dataIO.text.inputData1,h.dataIO.edit.inputData1,h.dataIO.button.inputData1] = sepia_construct_text_edit_button(panelParent,...
        'or Phase:',[],open_icon,pos,wratio);

    % row 3, col 1
    % input data2
    pos = [left(1) bottom(3) width height];
    [h.dataIO.text.inputData2,h.dataIO.edit.inputData2,h.dataIO.button.inputData2] = sepia_construct_text_edit_button(panelParent,...
        '    Magnitude:',[],open_icon,pos,wratio);

    % row 4, col 1
    % input data3
    pos = [left(1) bottom(4) width height];
    [h.dataIO.text.inputData3,h.dataIO.edit.inputData3,h.dataIO.button.inputData3] = sepia_construct_text_edit_button(panelParent,...
        '    Weights:',[],open_icon,pos,wratio);

    % row 5, col 1
    % input header
    pos = [left(1) bottom(5) width height];
    [h.dataIO.text.inputHeader,h.dataIO.edit.inputHeader,h.dataIO.button.inputHeader] = sepia_construct_text_edit_button(panelParent,...
        '    SEPIA header:',[],open_icon,pos,wratio);

    % row 4, col 2
    % invert pahse
    h.dataIO.checkbox.invertPhase = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','checkbox','String','Invert phase data',...
        'units','normalized','Position',[left(2) bottom(4) width height],...
        'backgroundcolor',get(h.fig,'color'));
    
    % row 5, col 2
    % brain extraction
    h.dataIO.checkbox.brainExtraction = uicontrol('Parent',h.StepsPanel.dataIO,...
        'Style','checkbox','String','FSL brain extraction (bet),',...
        'units','normalized','Position',[left(2) bottom(5) width*0.5 height],...
        'backgroundcolor',get(h.fig,'color'));
    
    % fractional threshold
    wratio = 0.5;
    [h.dataIO.text.fractionalThres,h.dataIO.edit.fractionalThres] = sepia_construct_text_edit(panelParent,...
        '-f',defaultFractionalThres,[left(2)+width*0.5 bottom(5) width*0.25 height],wratio);
    set(h.dataIO.text.fractionalThres, 'HorizontalAlignment','right');
    set(h.dataIO.edit.fractionalThres, 'Enable','off');
    
    % gradient threshold
    wratio = 0.5;
    [h.dataIO.text.gradientThres,h.dataIO.edit.gradientThres] = sepia_construct_text_edit(panelParent,...
        '-g',defaultGradientThres,[left(2)+width*0.75 bottom(5) width*0.25 height],wratio);
    set(h.dataIO.text.gradientThres, 'HorizontalAlignment','right');
    set(h.dataIO.edit.gradientThres, 'Enable','off');

    
%% set tooltips
set(h.dataIO.text.input,                'Tooltip',tooltip.dataIO.SEPIA.input);
set(h.dataIO.text.output,               'Tooltip',tooltip.dataIO.SEPIA.output);
set(h.dataIO.text.maskdir,              'Tooltip',tooltip.dataIO.SEPIA.maskdir);
set(h.dataIO.text.inputData1,           'Tooltip',tooltip.dataIO.SEPIA.inputData1);
set(h.dataIO.text.inputData2,           'Tooltip',tooltip.dataIO.SEPIA.inputData2);
set(h.dataIO.text.inputData3,           'Tooltip',tooltip.dataIO.SEPIA.inputData3);
set(h.dataIO.text.inputHeader,          'Tooltip',tooltip.dataIO.SEPIA.inputHeader);
set(h.dataIO.checkbox.invertPhase,      'Tooltip',tooltip.dataIO.SEPIA.invertPhase);
set(h.dataIO.checkbox.brainExtraction,	'Tooltip',tooltip.dataIO.SEPIA.BET);
set(h.dataIO.text.fractionalThres,      'Tooltip',tooltip.dataIO.SEPIA.fractThres);
set(h.dataIO.text.gradientThres,        'Tooltip',tooltip.dataIO.SEPIA.gradThres);

%% set callback function
% checkbox/edit pair
set(h.dataIO.checkbox.brainExtraction,  'Callback', {@CheckboxBrainExtraction_Callback,h});
% open/edit pari
set(h.dataIO.button.input,           	'Callback', {@ButtonOpen_Callback,h,'input'});
set(h.dataIO.button.output,         	'Callback', {@ButtonOpen_Callback,h,'output'});
set(h.dataIO.button.maskdir,          	'Callback', {@ButtonOpen_Callback,h,'mask'});
set(h.dataIO.button.inputData1,         'Callback', {@ButtonOpen_Callback,h,'inputdata1'});
set(h.dataIO.button.inputData2,         'Callback', {@ButtonOpen_Callback,h,'inputdata2'});
set(h.dataIO.button.inputData3,         'Callback', {@ButtonOpen_Callback,h,'inputdata3'});
set(h.dataIO.button.inputHeader,        'Callback', {@ButtonOpen_Callback,h,'header'});

set(h.dataIO.edit.fractionalThres,    	'Callback', {@EditInputMinMax_Callback,defaultFractionalThres,0,0,1});
set(h.dataIO.edit.gradientThres,        'Callback', {@EditInputMinMax_Callback,defaultGradientThres,0,-1,1});
    
end