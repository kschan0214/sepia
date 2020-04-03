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
% Date last modified:
%
%
function h = sepia_handle_panel_swi(hParent,h,position)
% set up method name displayed on GUI
methodName = {'SWI','SMWI'};

% Set parent of qsm panel
h.StepsPanel.swi = uipanel(hParent,...
    'Title','SWI','backgroundcolor',get(h.fig,'color'),...
    'position',[position(1) position(2) 0.95 0.35]);

%% design of this panel

    % text|popup pair: select method
    h.swi.text.swi = uicontrol('Parent',h.StepsPanel.swi,'Style','text','String','Method:',...
        'units','normalized','position',[0.01 0.85 0.15 0.1],...
        'HorizontalAlignment','left',...
        'backgroundcolor',get(h.fig,'color'),...
        'tooltip','Select SWI/SMWI');
    h.swi.popup.swi = uicontrol('Parent',h.StepsPanel.swi,'Style','popup',...
        'String',methodName,...
        'units','normalized','position',[0.31 0.85 0.4 0.1]) ; 
    
    % start button
    h.swi.button.start = uicontrol('Parent',h.StepsPanel.swi,...
        'Style','pushbutton','String','Start',...
        'units','normalized','position',[0.79 0.03 0.2 0.1],...
        'backgroundcolor','white');
    
% create control panel

% define position and size of all method panels
position_child = [0.01 0.15 0.95 0.6];

    % SWI    
    h = sepia_handle_panel_swi_SWI(h.StepsPanel.swi,h,position_child);

    % SMWI
    h = sepia_handle_panel_swi_SMWI(h.StepsPanel.swi,h,position_child);
    
    % in future, add panel of new method here

% set callback
set(h.swi.popup.swi,    'Callback', {@PopupSWI_Callback,h});
set(h.swi.button.start,	'Callback', {@PushbuttonStart_swi_Callback,h});
end

%% Callback function
% display corresponding method panel
function PopupSWI_Callback(source,eventdata,h)

% get selected QSM method
method = source.String{source.Value,1} ;

% switch off all panels
fields = fieldnames(h.swi.panel);
for kf = 1:length(fields)
    set(h.swi.panel.(fields{kf}),   'Visible','off');
end

% switch on target panel
switch method
    case 'SWI'
        set(h.swi.panel.SWI,        'Visible','on');

    case 'SMWI'
        set(h.swi.panel.SMWI,       'Visible','on');

    % in the future, add new method here

end

end

function PushbuttonStart_swi_Callback(source,eventdata,h)
% global h

% Disable the pushbutton to prevent double clicks
set(source,'Enable','off');

% get I/O GUI input
% option 2: input are NIfTI files
input(1).name = get(h.swi.dataIO.edit.inputData1, 	'String');
input(2).name = get(h.swi.dataIO.edit.inputData2,  	'String');

outputBasename  = get(h.swi.dataIO.edit.output,       	'String');

% get and create output directory
output_index = strfind(outputBasename, filesep);
outputDir = outputBasename(1:output_index(end));
% if the output directory does not exist then create the directory
if exist(outputDir,'dir') ~= 7
    mkdir(outputDir);
end

% create a new m file
% logFilename = [outputDir filesep 'sepia_log.m'];
% if exist(logFilename,'file') == 2
%     counter = 1;
%     while exist(logFilename,'file') == 2
%         suffix = ['_' num2str(counter)];
%         logFilename = [outputDir filesep 'sepia_log' suffix '.m'];
%         counter = counter + 1;
%     end
% end
% fid = fopen(logFilename,'w');
% create a new m file
configFilename = [outputDir filesep 'sepia_config.m'];
if exist(configFilename,'file') == 2
    counter = 1;
    while exist(configFilename,'file') == 2
        suffix = ['_' num2str(counter)];
        configFilename = [outputDir filesep 'sepia_config' suffix '.m'];
        counter = counter + 1;
    end
end
fid = fopen(configFilename,'w');

% input data
fprintf(fid,'input(1).name = ''%s'' ;\n',input(1).name);
fprintf(fid,'input(2).name = ''%s'' ;\n',input(2).name);

% output
fprintf(fid,'output_basename = ''%s'' ;\n',outputBasename);

% write algorithm parameters
switch h.swi.popup.swi.String{h.swi.popup.swi.Value,1}
    case 'SWI'
        fprintf(fid,'algorParam.swi.m = %s ;\n'                 ,get(h.swi.SWI.edit.m,  'String'));
        fprintf(fid,'algorParam.swi.threshold = %s ;\n'         ,get(h.swi.SWI.edit.threshold,  'String'));
        fprintf(fid,'algorParam.swi.filterSize = %s ;\n'        ,get(h.swi.SWI.edit.filterSize,  'String'));
        switch h.swi.SWI.popup.method.String{h.swi.SWI.popup.method.Value,1}
            case 'default'
                fprintf(fid,'algorParam.swi.method = ''%s'' ;\n'    ,'default');
            case 'multi-echo'
                fprintf(fid,'algorParam.swi.method = ''%s'' ;\n'    ,'multiecho');
        end
        fprintf(fid,'algorParam.swi.isPositive = %i ;\n'        ,get(h.swi.SWI.checkbox.positive,'Value'));
        fprintf(fid,'algorParam.swi.isNegative = %i ;\n'        ,get(h.swi.SWI.checkbox.negative,'Value'));
        fprintf(fid,'algorParam.swi.ismIP = %i ;\n'             ,get(h.swi.SWI.checkbox.mIP,'Value'));
        if get(h.swi.SWI.checkbox.mIP,'Value')
            fprintf(fid,'algorParam.swi.slice_mIP = %s ;\n' ,get(h.swi.SWI.edit.mIP,  'String'));
        end
        
        fprintf(fid,'\nSWIIOWrapper(input,output_basename,algorParam);\n');
        
    case 'SMWI'
        fprintf(fid,'algorParam.smwi.m = %s ;\n'                ,get(h.swi.SMWI.edit.m,  'String'));
        fprintf(fid,'algorParam.smwi.threshold = %s ;\n'        ,get(h.swi.SMWI.edit.threshold,  'String'));
        fprintf(fid,'algorParam.smwi.isParamagnetic = %i ;\n'   ,get(h.swi.SMWI.checkbox.paramagnetic,'Value'));
        fprintf(fid,'algorParam.smwi.isDiamagnetic = %i ;\n'    ,get(h.swi.SMWI.checkbox.diamagnetic,'Value'));
        fprintf(fid,'algorParam.smwi.ismIP = %i ;\n'            ,get(h.swi.SMWI.checkbox.mIP,'Value'));
        if get(h.swi.SMWI.checkbox.mIP,'Value')
            fprintf(fid,'algorParam.smwi.slice_mIP = %s ;\n' ,get(h.swi.SMWI.edit.mIP,  'String'));
        end
        
        fprintf(fid,'\nSMWIIOWrapper(input,output_basename,algorParam);\n');
end

fclose(fid);

% log command window display to a text file
logFilename = [outputDir filesep 'run_sepia.log'];
if exist(logFilename,'file') == 2
    counter = 1;
    while exist(logFilename,'file') == 2
        suffix = ['_' num2str(counter)];
        logFilename = [outputDir filesep 'run_sepia' suffix '.log'];
        counter = counter + 1;
    end
end
diary(logFilename)
    
try
    
    % run process
    run(logFilename);
    
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
    errorMessageFilename = [outputDir filesep 'run_sepia.error'];
    if exist(errorMessageFilename,'file') == 2
        counter = 1;
        while exist(errorMessageFilename,'file') == 2
            suffix = ['_' num2str(counter)];
            errorMessageFilename = [outputDir filesep 'run_sepia' suffix '.error'];
            counter = counter + 1;
        end
    end
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