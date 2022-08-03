%% h = sepia_handle_panel_swi(hParent,h,position)
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
% Description: This GUI function creates a panel for QSM method control
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 14 April 2019
% Date modified: 2 August 2022 (v1.1)
%
%
function h = sepia_handle_panel_swismwi(hParent,h,position)
% set up method name displayed on GUI
sepia_universal_variables;

%% Tooltips

% %% layout of the panel
% % define maximum level of options and spacing between options
% ncol    = 2; % 2 columns in the panel
% rspacing = 0.01;
% width   = (1-(ncol+1)*rspacing)/ncol;
% left    = (rspacing:width+rspacing:(width+rspacing)*ncol);

% % set up method name displayed on GUI
% methodName = {'SWI','SMWI'};

%% layout of the panel
% define maximum level of options and spacing between options
ncol    = 2; % 2 columns in the panel
rspacing = 0.01;
width   = (1-(ncol+1)*rspacing)/ncol;
left    = (rspacing:width+rspacing:(width+rspacing)*ncol);

%% Set parent of qsm panel
h.StepsPanel.swismwi = uipanel(hParent,...
    'Title','SWI/SMWI','backgroundcolor',get(h.fig,'color'),...
    'position',[position(1) position(2) 0.95 0.35]);

%% design of this panel

    height = 0.1;
    wratio = 0.5;

    % col 1
    % text|popup pair: select method
    [h.swismwi.text.method,h.swismwi.popup.method] = sepia_construct_text_popup(...
        h.StepsPanel.swismwi,'Method:', methodSWISMWIName, [left(1) 0.85 width height], wratio);
    
    % start button
    h.swismwi.button.start = uicontrol('Parent',h.StepsPanel.swismwi,...
        'Style','pushbutton','String','Start',...
        'units','normalized','position',[0.79 0.03 0.2 0.1],...
        'backgroundcolor','white');
    
% create control panel

% define position and size of all method panels
position_child = [0.01 0.15 0.95 0.6];

% construct all method panels
for k = 1:length(function_SWISMWI_method_panel)
    h = feval(function_SWISMWI_method_panel{k},h.StepsPanel.swismwi,h,position_child);
end

%     % SWI    
%     h = sepia_handle_panel_swi_SWI(h.StepsPanel.swi,h,position_child);
% 
%     % SMWI
%     h = sepia_handle_panel_swi_SMWI(h.StepsPanel.swi,h,position_child);
    
    % in future, add panel of new method here
    
%% set tooltip
set(h.swismwi.text.method,     'Tooltip','Select a SWI or SMWI method');

%% set callback
set(h.swismwi.popup.method,     'Callback', {@PopupSWI_Callback,h});
set(h.swismwi.button.start,     'Callback', {@PushbuttonStart_swi_Callback,h});
end

%% Callback function
% display corresponding method panel
function PopupSWI_Callback(source,eventdata,h)

sepia_universal_variables;

% get selected SWI/SMWI method
method = source.String{source.Value,1} ;

% switch off all panels first
fields = fieldnames(h.swismwi.panel);
for kf = 1:length(fields)
    set(h.swismwi.panel.(fields{kf}),    'Visible','off');
end

% switch on only target panel
for k = 1:length(methodSWISMWIName)
    if strcmpi(method,methodSWISMWIName{k})
        set(h.swismwi.panel.(fields{k}), 'Visible','on');
        break
    end
end

% % switch on target panel
% switch method
%     case 'SWI'
%         set(h.swi.panel.SWI,        'Visible','on');
% 
%     case 'SMWI'
%         set(h.swi.panel.SMWI,       'Visible','on');
% 
%     % in the future, add new method here
% 
% end

end

function PushbuttonStart_swi_Callback(source,eventdata,h)

% global h

sepia_universal_variables;

%%%%%%%%%%%% Step 1: preparation %%%%%%%%%%%%
% Disable the pushbutton to prevent double clicks
set(source,'Enable','off');

% get I/O GUI input
% option 2: input are NIfTI files
input(1).name = get(h.swismwi.dataIO.edit.inputData1, 	'String');
input(2).name = get(h.swismwi.dataIO.edit.inputData2,  	'String');
input(3).name = get(h.swismwi.dataIO.edit.inputHeader,	'String');

outputBasename  = get(h.swismwi.dataIO.edit.output,     'String');

% get and create output directory
output_index = strfind(outputBasename, filesep);
outputDir = outputBasename(1:output_index(end));
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

% use current time as unique identifier
identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');

% create a new m file
configFilename = fullfile(outputDir, ['sepia_config_' identifier '.m']);
while exist(configFilename,'file') == 2
    % update current time as unique identifier
    identifier = datestr(datetime('now'),'yymmddHHMMSSFFF');
    configFilename = fullfile(outputDir, ['sepia_config_' identifier '.m']);
end
% configFilename = [outputDir filesep 'sepia_config.m'];
% if exist(configFilename,'file') == 2
%     counter = 1;
%     while exist(configFilename,'file') == 2
%         suffix = ['_' num2str(counter)];
%         configFilename = [outputDir filesep 'sepia_config' suffix '.m'];
%         counter = counter + 1;
%     end
% end
fid = fopen(configFilename,'w');

%%%%%%%%%%%% Step 2: Write config file %%%%%%%%%%%% 
% specify config file version
fprintf(fid,'%% This file is generated by SEPIA version: %s\n',SEPIA_version);
% general path
fprintf(fid,'%% add general Path\n');
fprintf(fid,'sepia_addpath;\n\n');

fprintf(fid,'%% Input/Output filenames\n');
% input data
fprintf(fid,'input(1).name = ''%s'' ;\n',input(1).name);
fprintf(fid,'input(2).name = ''%s'' ;\n',input(2).name);
fprintf(fid,'input(3).name = ''%s'' ;\n',input(3).name);

% output
fprintf(fid,'output_basename = ''%s'' ;\n',outputBasename);

fprintf(fid,'%% SWI/SMWI algorithm parameters\n');
% set parameters for selected method
print_method_popup_and_eval(fid, '.swismwi.method', h.swismwi.popup.method, methodSWISMWIName, config_SWISMWI_function, h);

% Determine application based on Tab
fprintf(fid,'\nSWISMWIIOWrapper(input,output_basename,algorParam);\n');

% % write algorithm parameters
% switch h.swismwi.popup.swi.String{h.swismwi.popup.swi.Value,1}
%     case 'SWI'
%         fprintf(fid,'algorParam.swi.m = %s ;\n'                 ,get(h.swi.SWI.edit.m,  'String'));
%         fprintf(fid,'algorParam.swi.threshold = %s ;\n'         ,get(h.swi.SWI.edit.threshold,  'String'));
%         fprintf(fid,'algorParam.swi.filterSize = %s ;\n'        ,get(h.swi.SWI.edit.filterSize,  'String'));
%         switch h.swi.SWI.popup.method.String{h.swi.SWI.popup.method.Value,1}
%             case 'default'
%                 fprintf(fid,'algorParam.swi.method = ''%s'' ;\n'    ,'default');
%             case 'multi-echo (testing)'
%                 fprintf(fid,'algorParam.swi.method = ''%s'' ;\n'    ,'multiecho (testing)');
%         end
%         fprintf(fid,'algorParam.swi.isPositive = %i ;\n'        ,get(h.swi.SWI.checkbox.positive,'Value'));
%         fprintf(fid,'algorParam.swi.isNegative = %i ;\n'        ,get(h.swi.SWI.checkbox.negative,'Value'));
%         fprintf(fid,'algorParam.swi.ismIP = %i ;\n'             ,get(h.swi.SWI.checkbox.mIP,'Value'));
%         if get(h.swi.SWI.checkbox.mIP,'Value')
%             fprintf(fid,'algorParam.swi.slice_mIP = %s ;\n' ,get(h.swi.SWI.edit.mIP,  'String'));
%         end
%         
%         fprintf(fid,'\nSWIIOWrapper(input,output_basename,algorParam);\n');
%         
%     case 'SMWI'
%         fprintf(fid,'algorParam.smwi.m = %s ;\n'                ,get(h.swi.SMWI.edit.m,  'String'));
%         fprintf(fid,'algorParam.smwi.threshold = %s ;\n'        ,get(h.swi.SMWI.edit.threshold,  'String'));
%         fprintf(fid,'algorParam.smwi.isParamagnetic = %i ;\n'   ,get(h.swi.SMWI.checkbox.paramagnetic,'Value'));
%         fprintf(fid,'algorParam.smwi.isDiamagnetic = %i ;\n'    ,get(h.swi.SMWI.checkbox.diamagnetic,'Value'));
%         fprintf(fid,'algorParam.smwi.ismIP = %i ;\n'            ,get(h.swi.SMWI.checkbox.mIP,'Value'));
%         if get(h.swi.SMWI.checkbox.mIP,'Value')
%             fprintf(fid,'algorParam.smwi.slice_mIP = %s ;\n' ,get(h.swi.SMWI.edit.mIP,  'String'));
%         end
%         
%         fprintf(fid,'\nSMWIIOWrapper(input,output_basename,algorParam);\n');
% end

fclose(fid);

% log command window display to a text file
% logFilename = fullfile(outputDir, 'run_sepia.log'); 
logFilename = fullfile(outputDir, ['run_sepia.log' identifier]);
% if exist(logFilename,'file') == 2
%     counter = 1;
%     while exist(logFilename,'file') == 2
%         suffix = ['_' num2str(counter)];
%         logFilename = [outputDir filesep 'run_sepia' suffix '.log'];
%         counter = counter + 1;
%     end
% end
diary(logFilename)
    
try
    
    % run process
    run(configFilename);
    
    % turn off the log
    diary off

    % re-enable the pushbutton
    set(source,'Enable','on');

catch ME
    % re-enable the start button before displaying the error
    set(source,'Enable','on');
    
     % close log file
    disp('There was an error! Please check the command window/error message file for more information.');
    diary off
    
    % open a new text file for error message
    errorMessageFilename = fullfile(outputDir, ['run_sepia.error' identifier]);
%     errorMessageFilename = [outputDir filesep 'run_sepia.error'];
%     if exist(errorMessageFilename,'file') == 2
%         counter = 1;
%         while exist(errorMessageFilename,'file') == 2
%             suffix = ['_' num2str(counter)];
%             errorMessageFilename = [outputDir filesep 'run_sepia' suffix '.error'];
%             counter = counter + 1;
%         end
%     end
    fid = fopen(errorMessageFilename,'w');
    fprintf(fid,'The identifier was:\n%s\n\n',ME.identifier);
    fprintf(fid,'The message was:\n\n');
    msgString = getReport(ME,'extended','hyperlinks','off');
    fprintf(fid,'%s',msgString);
    fclose(fid);
    
    % rethrow the error message to command window
    rethrow(ME);
end

        

end

function print_method_popup_and_eval(fid, str_pattern, action_handle, popup_list, config_list, h)

% set parameters for selected method
method = sepia_print_popup_as_string(fid,str_pattern,action_handle);
for k = 1:length(popup_list)
    if strcmpi(method,popup_list{k})
        feval(config_list{k},h,'set',fid);
    end
end
    
end