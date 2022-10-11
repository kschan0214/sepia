%% h = sepia_handle_panel_utility_get_header(hParent,hFig,h,position)
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
% Description: This GUI function creates a panel for the utility function
% 'Get qsm_hub header'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 24 May 2018
% Date modified: 12 June 2018
% Date modified: 24 May 2019
%
%
function h = sepia_handle_panel_analysis_segmentation_AHEAD(hParent,h,position)

open_icon = imread('folder@0,3x.jpg');
open_icon = imresize(open_icon,[1 1]*16);

%% layout of the panel
nrow        = 20;
rspacing    = 0.02;
ncol        = 1;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);

%%
% set Parent of all related controls
h.Analysis.panel.Segmentation_AHEAD = uipanel(hParent,'Title','Segmentation - AHEAD atlas',...
    'Position',position,...
    'backgroundcolor',get(h.fig,'color'),'Visible','off');
    
    panelParent = h.Analysis.panel.Segmentation_AHEAD;
    
    wratio = [0.3,0.65,0.05];
    % Option 1: NIfTI file input
    h.Analysis.segmentation.AHEAD.text.Option1 = uicontrol('Parent',panelParent ,'Style','text','units','normalized', 'HorizontalAlignment','left', 'backgroundcolor',get(gcf,'color'),'FontWeight','bold',...
        'String','Input Option 1','position',[left(1) bottom(1) width height]);
    
    pos = [left(1) bottom(2) width height];
    [h.Analysis.segmentation.AHEAD.text.greInput,...
     h.Analysis.segmentation.AHEAD.edit.greInput,...
     h.Analysis.segmentation.AHEAD.button.greInput] = sepia_construct_text_edit_button(panelParent,...
        'Select a 3D GRE magnitude NIfTI file:',[],open_icon,pos,wratio);
    
    pos = [left(1) bottom(3) width height];
    [h.Analysis.segmentation.AHEAD.text.greMaskInput,...
     h.Analysis.segmentation.AHEAD.edit.greMaskInput,...
     h.Analysis.segmentation.AHEAD.button.greMaskInput] = sepia_construct_text_edit_button(panelParent,...
        'Select a GRE mask NIfTI file:',[],open_icon,pos,wratio);
    
    pos = [left(1) bottom(4) width height];
    [h.Analysis.segmentation.AHEAD.text.t1wInput,...
     h.Analysis.segmentation.AHEAD.edit.t1wInput,...
     h.Analysis.segmentation.AHEAD.button.t1wInput] = sepia_construct_text_edit_button(panelParent,...
        'Select a T1w NIfTI file:',[],open_icon,pos,wratio);
    
    pos = [left(1) bottom(5) width height];
    [h.Analysis.segmentation.AHEAD.text.t1wMaskInput,...
     h.Analysis.segmentation.AHEAD.edit.t1wMaskInput,...
     h.Analysis.segmentation.AHEAD.button.t1wMaskInput] = sepia_construct_text_edit_button(panelParent,...
        'Select a T1w mask NIfTI file:',[],open_icon,pos,wratio);
    
    pos = [left(1) bottom(6) width height];
    [h.Analysis.segmentation.AHEAD.text.chiInput,...
     h.Analysis.segmentation.AHEAD.edit.chiInput,...
     h.Analysis.segmentation.AHEAD.button.chiInput] = sepia_construct_text_edit_button(panelParent,...
        'Select a Chimap NIfTI file:',[],open_icon,pos,wratio);
    
    % Option 2: Transformation input
    h.Analysis.segmentation.AHEAD.text.Option2 = uicontrol('Parent',panelParent ,'Style','text','units','normalized', 'HorizontalAlignment','left', 'backgroundcolor',get(gcf,'color'),'FontWeight','bold',...
        'String','Input Option 2','position',[left(1) bottom(7) width height]);
    pos = [left(1) bottom(8) width height];
    [h.Analysis.segmentation.AHEAD.text.greInput2,...
     h.Analysis.segmentation.AHEAD.edit.greInput2,...
     h.Analysis.segmentation.AHEAD.button.greInput2] = sepia_construct_text_edit_button(panelParent,...
        'Select a GRE magnitude NIfTI file:',[],open_icon,pos,wratio);
    pos = [left(1) bottom(9) width height];
    [h.Analysis.segmentation.AHEAD.text.gre2T1wMat,...
     h.Analysis.segmentation.AHEAD.edit.gre2T1wMat,...
     h.Analysis.segmentation.AHEAD.button.gre2T1wMat] = sepia_construct_text_edit_button(panelParent,...
        'Select a GRE-to-T1w rigid-body transformation:',[],open_icon,pos,wratio);
    pos = [left(1) bottom(10) width height];
    [h.Analysis.segmentation.AHEAD.text.t1w2TemplateMat,...
     h.Analysis.segmentation.AHEAD.edit.t1w2TemplateMat,...
     h.Analysis.segmentation.AHEAD.button.t1w2TemplateMat] = sepia_construct_text_edit_button(panelParent,...
        'Select a T1w-to-Atlas affine transformation:',[],open_icon,pos,wratio);
    pos = [left(1) bottom(11) width height];
    [h.Analysis.segmentation.AHEAD.text.t1w2TemplateNii,...
     h.Analysis.segmentation.AHEAD.edit.t1w2TemplateNii,...
     h.Analysis.segmentation.AHEAD.button.t1w2TemplateNii] = sepia_construct_text_edit_button(panelParent,...
        'Select a T1w-to-Atlas Inverse Wrap NIfTI file:',[],open_icon,pos,wratio);
    
    % output
    pos = [left(1) bottom(13) width height];
    [h.Analysis.segmentation.AHEAD.text.outputDir,...
     h.Analysis.segmentation.AHEAD.edit.outputDir,...
     h.Analysis.segmentation.AHEAD.button.outputDir] = sepia_construct_text_edit_button(panelParent,...
        'Ourput directory:',pwd,open_icon,pos,wratio);
    % correct bias field option
    pos = [left(1) bottom(14) width height];
    h.Analysis.segmentation.AHEAD.checkbox.biasCorr = uicontrol('Parent',panelParent,'backgroundcolor',get(h.fig,'color'),'Style','checkbox','units','normalized',...
        'String','Correct bias field on input images','Position',pos);
  
    % run
    pos = [0.79 bottom(end) 0.2 height*2];
    h.Analysis.segmentation.AHEAD.button.start = uicontrol('Parent',panelParent,'Style','pushbutton','backgroundcolor','white','units','normalized',...
        'String','Start', 'position',pos);
%     
%% set callback functions
set(h.Analysis.segmentation.AHEAD.button.greInput,           'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'nitfi',h.Analysis.segmentation.AHEAD.edit.greInput,[1 2]});
set(h.Analysis.segmentation.AHEAD.button.greMaskInput,       'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'nitfi',h.Analysis.segmentation.AHEAD.edit.greMaskInput,[1 2]});
set(h.Analysis.segmentation.AHEAD.button.t1wInput,           'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'nitfi',h.Analysis.segmentation.AHEAD.edit.t1wInput,[1 2]});
set(h.Analysis.segmentation.AHEAD.button.t1wMaskInput,       'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'nitfi',h.Analysis.segmentation.AHEAD.edit.t1wMaskInput,[1 2]});
set(h.Analysis.segmentation.AHEAD.button.chiInput,           'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'nitfi',h.Analysis.segmentation.AHEAD.edit.chiInput,[1 2]});

set(h.Analysis.segmentation.AHEAD.button.greInput2,        	'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'nitfi',h.Analysis.segmentation.AHEAD.edit.greInput2,[1 1]});
set(h.Analysis.segmentation.AHEAD.button.gre2T1wMat,         'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'mat',h.Analysis.segmentation.AHEAD.edit.gre2T1wMat,[1 1]});
set(h.Analysis.segmentation.AHEAD.button.t1w2TemplateMat,	'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'mat',h.Analysis.segmentation.AHEAD.edit.t1w2TemplateMat,[1 1]});
set(h.Analysis.segmentation.AHEAD.button.t1w2TemplateNii,	'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'nitfi',h.Analysis.segmentation.AHEAD.edit.t1w2TemplateNii,[1 1]});
set(h.Analysis.segmentation.AHEAD.button.outputDir,          'Callback', {@ButtonOpen_Analysis_segmentation_Callback,h,'dir',h.Analysis.segmentation.AHEAD.edit.outputDir,[0]});

set(h.Analysis.segmentation.AHEAD.button.start,             'Callback', {@PushbuttonStart_Analysis_segmentation_Callback,h});
end

%% Callback functions
% 'open' button callback
function ButtonOpen_Analysis_segmentation_Callback(source,eventdata,h,field,field_handle,clearField)
% get input file/directory for getHeader utility function
% global h

switch field
    case 'nitfi'
        % only read NIfTI file for mask
        [niftiFileName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select a NIFTI file');
        
        % set edit field
        if pathDir ~= 0
            set(field_handle,'String',fullfile(pathDir,niftiFileName));
        end
        
    case 'mat'
        % only read NIfTI file for mask
        [matName,pathDir] = uigetfile({'*.mat'},'Select an ANTs Tranformation file');
        
        % set edit field
        if pathDir ~= 0
            set(field_handle,'String',fullfile(pathDir,matName));
        end
        
    case 'dir'
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(field_handle,'String', fullfile(pathDir));
        end
        
end

if clearField(1)
    switch clearField(2)
        case 1
            set(h.Analysis.segmentation.AHEAD.edit.greInput,'String','');
            set(h.Analysis.segmentation.AHEAD.edit.greMaskInput,'String','');
            set(h.Analysis.segmentation.AHEAD.edit.t1wInput,'String','');
            set(h.Analysis.segmentation.AHEAD.edit.t1wMaskInput,'String','');
            set(h.Analysis.segmentation.AHEAD.edit.chiInput,'String','');
            
        case 2
            set(h.Analysis.segmentation.AHEAD.edit.greInput2,'String','');
            set(h.Analysis.segmentation.AHEAD.edit.gre2T1wMat,'String','');
            set(h.Analysis.segmentation.AHEAD.edit.t1w2TemplateMat,'String','');
            set(h.Analysis.segmentation.AHEAD.edit.t1w2TemplateNii,'String','');
    end
end

end

function PushbuttonStart_Analysis_segmentation_Callback(source,eventdata,h)
% Callback function to detect and save header functino for sepia
sepia_universal_variables;
% global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

% Display message
disp('Perform image registration to get atlas labels...');

% get input from GUI
% option 1
input.gre       = get(h.Analysis.segmentation.AHEAD.edit.greInput,       'String');
input.greMask 	= get(h.Analysis.segmentation.AHEAD.edit.greMaskInput, 	'String');
input.t1w       = get(h.Analysis.segmentation.AHEAD.edit.t1wInput,       'String');
input.t1wMask  	= get(h.Analysis.segmentation.AHEAD.edit.t1wMaskInput, 	'String');
input.chi       = get(h.Analysis.segmentation.AHEAD.edit.chiInput,       'String');
if ~isempty(input.greMask)
    isPathway1 = true;
else
    % option 2  
    input.gre               = get(h.Analysis.segmentation.AHEAD.edit.greInput2,      'String');
    input.gre2T1wMat        = get(h.Analysis.segmentation.AHEAD.edit.gre2T1wMat, 	'String');
    input.t1w2TemplateMat 	= get(h.Analysis.segmentation.AHEAD.edit.t1w2TemplateMat,'String');
    input.t1w2TemplateiWrap	= get(h.Analysis.segmentation.AHEAD.edit.t1w2TemplateNii,'String');
    isPathway1 = false;
end
% output
outputDir       = get(h.Analysis.segmentation.AHEAD.edit.outputDir,       	'String');
% if the output directory does not exist then create the directory
if isempty(outputDir); outputDir = fullfile(pwd,'segmentation_sepia');end
if exist(outputDir,'dir') ~= 7; mkdir(outputDir); end

% create a new m file
identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
configFilename = [outputDir filesep 'sepia_segmentation_config' identifier '.m'];
if exist(configFilename,'file') == 2
    while exist(configFilename,'file') == 2
        % update current time as unique identifier
        identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
        configFilename = [outputDir filesep 'sepia_segmentation_config' identifier '.m'];
    end
end
fid = fopen(configFilename,'w');

%%%%%%%%%%%% Step 2: Write config file %%%%%%%%%%%% 
% specify config file version
fprintf(fid,'%% This file is generated by SEPIA version: %s\n',SEPIA_version);
% general path
fprintf(fid,'%% add general Path\n');
fprintf(fid,'sepia_addpath;\n\n');

fprintf(fid,'%% Input/Output filenames\n');

fprintf(fid,'input = struct();\n');
% input data
if isPathway1
    fprintf(fid,'input.gre      = ''%s'' ;\n',	input.gre);
    fprintf(fid,'input.greMask  = ''%s'' ;\n', 	input.greMask);
    fprintf(fid,'input.t1w      = ''%s'' ;\n', 	input.t1w);
    fprintf(fid,'input.t1wMask  = ''%s'' ;\n', 	input.t1wMask);
    fprintf(fid,'input.chi      = ''%s'' ;\n', 	input.chi);
else
    fprintf(fid,'input.gre                  = ''%s'' ;\n',	input.gre);
    fprintf(fid,'input.gre2T1wMat           = ''%s'' ;\n',	input.gre2T1wMat);
    fprintf(fid,'input.t1w2TemplateMat      = ''%s'' ;\n',	input.t1w2TemplateMat);
    fprintf(fid,'input.t1w2TemplateiWrap    = ''%s'' ;\n',	input.t1w2TemplateiWrap);
end
% output
fprintf(fid,'output_dir     = ''%s'' ;\n',outputDir);
fprintf(fid,'%% General algorithm parameters\n');
fprintf(fid,'algorParam = struct();\n');

% inform input method
if isPathway1; fprintf(fid,'algorParam.mode = 1;\n'); else; fprintf(fid,'algorParam.mode = 2;\n'); end
    
% is bias corr
sepia_print_checkbox_value(fid,'.isBiasFieldCorr',h.Analysis.segmentation.AHEAD.checkbox.biasCorr);

% Determine application based on Tab
fprintf(fid,'\nget_AHEAD_atlas_labels(input,output_dir,algorParam);\n');
fclose(fid);

try 
    

    logFilename = fullfile(outputDir, ['run_sepia_segmentation.log' identifier]);
    diary(logFilename)
    
    % run process
    run(configFilename);
    % re-enable the pushbutton
    set(source,'Enable','on');
    
    % turn off the log
    diary off
    
catch ME
    % close log file
    disp('There was an error! Please check the command window/error message file for more information.');
    diary off
    
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end

% re-enable the pushbutton 
set(source,'Enable','on');

end