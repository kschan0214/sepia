%% h = sepia_handle_panel_utility_convert_realImaginary2phase(hParent,h,position)
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
% 'Convert GE real/imaginary images to phase images'
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 22 September 2022
% Date modified:
%
%
function h = sepia_handle_panel_utility_convert_realImaginary2phase(hParent,h,position)

open_icon = imread('folder@0,3x.jpg');
open_icon = imresize(open_icon,[1 1]*16);

%% layout of the panel
nrow        = 6;
rspacing    = 0.02;
ncol        = 1;
cspacing    = 0.01;
[height,bottom,width,left] = sepia_layout_measurement(nrow,rspacing,ncol,cspacing);
    
wratio = [0.3,0.65,0.05];

% set Parent of all related controls
h.Utility.panel.realimag2phase = uipanel(hParent,'Title','Convert GE real/imaginary images to phase images',...
    'Position',position,...
    'backgroundcolor',get(h.fig,'color'),...
    'visible','off');
    
    panelParent = h.Utility.panel.realimag2phase;

    % BIDS file input
    pos = [left(1) bottom(1) width height];
    [h.Utility.realimag2phase.text.bidsInput,h.Utility.realimag2phase.edit.bidsInput,h.Utility.realimag2phase.button.bidsInput] = sepia_construct_text_edit_button(panelParent,...
        '(Option 1) BIDS directory:',[],open_icon,pos,wratio);

    % real file input
    pos = [left(1) bottom(2) width height];
    [h.Utility.realimag2phase.text.realInput,h.Utility.realimag2phase.edit.realInput,h.Utility.realimag2phase.button.realInput] = sepia_construct_text_edit_button(panelParent,...
        'or (Option 2) Real NIfTI image:',[],open_icon,pos,wratio);
    
    % imaginary file
    pos = [left(1) bottom(3) width height];
    [h.Utility.realimag2phase.text.imagInput,h.Utility.realimag2phase.edit.imagInput,h.Utility.realimag2phase.button.imagInput] = sepia_construct_text_edit_button(panelParent,...
        'and Imaginary NIfTI image:',[],open_icon,pos,wratio);
    
    pos = [left(1) bottom(4) width height];
    h.Utility.realimag2phase.checkbox.isCorrectInterslice = uicontrol('Parent',h.Utility.panel.realimag2phase ,...
        'Style','checkbox','String','Correct interslice phase polarity',...
        'units','normalized','Position',pos,...
        'backgroundcolor',get(h.fig,'color'),'Value',true);
    
    % output directory
    pos = [left(1) bottom(5) width height];
    [h.Utility.realimag2phase.text.outputDir,h.Utility.realimag2phase.edit.outputDir,h.Utility.realimag2phase.button.outputDir] = sepia_construct_text_edit_button(panelParent,...
        'Output directory with filename prefix',[],open_icon,pos,wratio);
    
    % run
    pos = [0.89 bottom(6) 0.1 height];
    h.Utility.realimag2phase.button.run = uicontrol('Parent',h.Utility.panel.realimag2phase,...
        'Style','pushbutton','String','Run',...
        'units','normalized','position',pos,...
        'backgroundcolor','white');

%% set callbacks
set(h.Utility.realimag2phase.button.imagInput,  	'Callback',     	{@ButtonOpen_Utility_realimag2phase_Callback,h,'imag'});
set(h.Utility.realimag2phase.button.realInput,      'Callback',         	{@ButtonOpen_Utility_realimag2phase_Callback,h,'real'});
set(h.Utility.realimag2phase.button.bidsInput,      'Callback',        	{@ButtonOpen_Utility_realimag2phase_Callback,h,'bidsdir'});
set(h.Utility.realimag2phase.button.outputDir,      'Callback',        	{@ButtonOpen_Utility_realimag2phase_Callback,h,'output'});
set(h.Utility.realimag2phase.button.run,            'Callback',       	{@PushbuttonRun_Utility_realimag2phase_Callback,h});
    
end

%% Callback functions
function ButtonOpen_Utility_realimag2phase_Callback(source,eventdata,h,field)
% get input file/directory for getHeader utility function
% global h
prefix = 'Sepia';
switch field
    case 'imag'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select NIfTI file');

        if pathDir ~= 0
            set(h.Utility.realimag2phase.edit.imagInput,    'String',fullfile(pathDir,nitfiName));
            % automatically set default output field
%             set(h.Utility.realimag2phase.edit.outputDir,     'String',[pathDir filesep 'output']);
        end

    case 'real'
        % read NIfTI file 
        [nitfiName,pathDir] = uigetfile({'*.nii;*.nii.gz','NIfTI file (*.nii,*.nii.gz)'},'Select NIfTI file');

        if pathDir ~= 0
            set(h.Utility.realimag2phase.edit.realInput,    'String',fullfile(pathDir,nitfiName));
            % automatically set default output field
            set(h.Utility.realimag2phase.edit.outputDir,  	'String', fullfile(pathDir,'output',prefix));
            set(h.Utility.realimag2phase.edit.bidsInput,    'String',[]);
        end
        
    case 'bidsdir'
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.Utility.realimag2phase.edit.bidsInput,     'String',pathDir);
            
            % automatically set default output field
            set(h.Utility.realimag2phase.edit.outputDir,  	'String', fullfile(pathDir,'output',prefix));
            set(h.Utility.realimag2phase.edit.realInput,    'String',[]);
            set(h.Utility.realimag2phase.edit.imagInput,    'String',[]);
        end
        
    case 'output'
        
        % get directory for output
        pathDir = uigetdir;

        if pathDir ~= 0
            set(h.Utility.realimag2phase.edit.outputDir,     'String',fullfile(pathDir,prefix));
        end
end

end

function PushbuttonRun_Utility_realimag2phase_Callback(source,event,h)
% Callback function to get lateral ventricle mask

% global h

% Disable the pushbutton to prevent doubel click
set(source,'Enable','off');

bidsDir                     = get(h.Utility.realimag2phase.edit.bidsInput, 'String');
real_img_fn                 = get(h.Utility.realimag2phase.edit.realInput, 'String');
imag_img_fn                 = get(h.Utility.realimag2phase.edit.imagInput, 'String');
isCorrInterSlicePolarity    = get(h.Utility.realimag2phase.checkbox.isCorrectInterslice, 'Value');
outputprefix            	= get(h.Utility.realimag2phase.edit.outputDir, 'String');

try
    disp('Converting real/imaginary images to phase images');
    
    if isempty(bidsDir) && and(~isempty(real_img_fn),~isempty(imag_img_fn))
        % NIfTI input
        img_real = load_nii_img_only(real_img_fn);
        img_imag = load_nii_img_only(imag_img_fn);

        if isCorrInterSlicePolarity
            img_real(:,:,2:2:end,:) = -img_real(:,:,2:2:end,:);
            img_imag(:,:,2:2:end,:) = -img_imag(:,:,2:2:end,:);
        end
        img_phase = angle(complex(img_real, img_imag));
        
        output_fn = strcat(outputprefix,'_part-phase.nii.gz');
        save_nii_img_only(real_img_fn,output_fn, img_phase);
        
    elseif ~isempty(bidsDir)
        % directory input
        real_list       = dir(fullfile(bidsDir,'*part-real*.nii*'));
        imaginary_list  = dir(fullfile(bidsDir,'*part-imaginary*.nii*'));
        
        % multi-echo 
        if numel(real_list) > 1
            for ke = 1:numel(real_list)
                img_real = load_nii_img_only(fullfile(real_list(ke).folder,fullfile(real_list(ke).name)));
                img_imag = load_nii_img_only(fullfile(imaginary_list(ke).folder,fullfile(imaginary_list(ke).name)));
                
                if isCorrInterSlicePolarity
                    img_real(:,:,2:2:end,:) = -img_real(:,:,2:2:end,:);
                    img_imag(:,:,2:2:end,:) = -img_imag(:,:,2:2:end,:);
                end
                img_phase = angle(complex(img_real, img_imag));
                
                output_fn = strcat(outputprefix,'_part-phase_echo-',num2str(ke),'.nii.gz');
                save_nii_img_only(fullfile(real_list(ke).folder,fullfile(real_list(ke).name)),output_fn, img_phase);
            end
        end
        
    end
    
    disp('Done!')

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    error(ME.message);
end
    

% Disable the pushbutton to prevent doubel click
set(source,'Enable','on');

end