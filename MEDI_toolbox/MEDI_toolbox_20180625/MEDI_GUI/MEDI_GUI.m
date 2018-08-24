function varargout = MEDI_GUI(varargin)
% MEDI_GUI MATLAB code for MEDI_GUI.fig
%      MEDI_GUI, by itself, creates a new MEDI_GUI or raises the existing
%      singleton*.
%
%      H = MEDI_GUI returns the handle to a new MEDI_GUI or the handle to
%      the existing singleton*.
%
%      MEDI_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEDI_GUI.M with the given input arguments.
%
%      MEDI_GUI('Property','Value',...) creates a new MEDI_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEDI_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEDI_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEDI_GUI

% Last Modified by GUIDE v2.5 02-Nov-2017 13:28:35

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MEDI_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MEDI_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MEDI_GUI is made visible.
function MEDI_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEDI_GUI (see VARARGIN)

% Choose default command line output for MEDI_GUI
handles.output = hObject;
handles.unwrapSelection = '';
handles.varin = '';
handles.loadedDataFolder = '';
handles.filedir = '';
handles.resultFolder = '';
handles.SMVEnable = 0;
handles.currPath = '';
set(handles.uipanel13,'selectedobject',handles.advRB);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEDI_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MEDI_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in BrowseDataButton.
function BrowseDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ExPath= uigetdir();

% ExPath = fullfile(FilePath, FileName);
set(handles.DataLocationEditable,'string',ExPath);


function DataLocationEditable_Callback(hObject, eventdata, handles)
% hObject    handle to DataLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DataLocationEditable as text
%        str2double(get(hObject,'String')) returns contents of DataLocationEditable as a double


% --- Executes during object creation, after setting all properties.
function DataLocationEditable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DataLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over BrowseDataButton.
function BrowseDataButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to BrowseDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Loadbutton.
function Loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FilePath = get(handles.DataLocationEditable,'string');
set(handles.statusText,'String','Busy Loading Data...  Please Wait');
drawnow
[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir] = Read_DICOM(FilePath);
set(handles.statusText,'String','Done Loading Data');
handles.loadedDataFolder = FilePath;
assignin('base','iField',iField);
assignin('base','voxel_size',voxel_size);
assignin('base','matrix_size',matrix_size);
assignin('base','CF',CF);
assignin('base','delta_TE',delta_TE);
assignin('base','TE',TE);
assignin('base','B0_dir',B0_dir);
guidata(hObject, handles);


% --- Executes on button press in FitButton.
function FitButton_Callback(hObject, eventdata, handles)
% hObject    handle to FitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iField = evalin('base','iField');
TE = evalin('base','TE');
set(handles.statusText,'String','Busy Fitting...');
drawnow
[iFreq_raw N_std] = Fit_ppm_complex_TE(iField,TE);
set(handles.statusText,'String','Fitting Done');
assignin('base','iFreq_raw',iFreq_raw);
assignin('base','N_std',N_std);

% % --- Executes on button press in pushbutton5.
% function pushbutton5_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton5 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles.Mask = genMask(handles.iField, handles.voxel_size);


% --- Executes on button press in BrowseMaskButton.
function BrowseMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileNameMask,FilePathMask]= uigetfile();
MaskFile = fullfile(FilePathMask, FileNameMask);
set(handles.MaskLocationEditable,'string',MaskFile);


function MaskLocationEditable_Callback(hObject, eventdata, handles)
% hObject    handle to MaskLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaskLocationEditable as text
%        str2double(get(hObject,'String')) returns contents of MaskLocationEditable as a double


% --- Executes during object creation, after setting all properties.
function MaskLocationEditable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaskLocationEditable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'Auto'
        iField = evalin('base','iField');
        voxel_size = evalin('base','voxel_size');
        mask = get(eventdata.NewValue,'Tag');
        Mask = genMask(iField, voxel_size);
        assignin('base','mask',mask);
        assignin('base','Mask',Mask);
        set(handles.MaskLocationEditable,'Enable','off');
        set(handles.BrowseMaskButton,'Enable','off');
        set(handles.LoadMaskButton,'Enable','off');
        set(handles.statusText,'String','Mask Automatically Generated');
        drawnow
    case 'BET'
        iField = evalin('base','iField');
        voxel_size = evalin('base','voxel_size');
        matrix_size = evalin('base','matrix_size');
        mask = get(eventdata.NewValue,'Tag');
        TE = evalin('base','TE');
        try
            iMag = evalin('base','iMag');
        catch
            iMag = sqrt(sum(abs(iField).^2,4));
            assignin('base','iMag',iMag);
        end
        Mask = BET(iMag, matrix_size, voxel_size);
        R2s = arlo(TE, abs(iField));
        Mask_CSF = extract_CSF(R2s, Mask, voxel_size);
        assignin('base','mask',mask);
        assignin('base','Mask',Mask);
        assignin('base','R2s',R2s);
        assignin('base','Mask_CSF',Mask_CSF);
        set(handles.MaskLocationEditable,'Enable','off');
        set(handles.BrowseMaskButton,'Enable','off');
        set(handles.LoadMaskButton,'Enable','off');
        set(handles.statusText,'String','Mask Automatically Generated');
        drawnow
    case 'UserSelect'
        mask = get(eventdata.NewValue,'Tag');
        assignin('base','mask',mask);
        set(handles.MaskLocationEditable,'Enable','on');
        set(handles.BrowseMaskButton,'Enable','on');
        set(handles.LoadMaskButton,'Enable','on');
    otherwise
        disp(':(');
end


% --- Executes on button press in UnwrapButton.
function UnwrapButton_Callback(hObject, eventdata, handles)
% Laplacian
% hObject    handle to UnwrapButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
switch get(get(handles.uipanel9,'SelectedObject'),'Tag')
    case'LaplacianRadioButton'
        iFreq_raw = evalin('base','iFreq_raw');
        matrix_size = evalin('base','matrix_size');
        voxel_size = evalin('base','voxel_size');
        set(handles.statusText,'String','Busy Performing Laplacian Unwrapping');
        drawnow
        iFreq = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size);
        set(handles.statusText,'String','Done Performing Laplacian Unwrapping');
        assignin('base','iFreq',iFreq);
        
    case'RGRadioButton'
        iField = evalin('base','iField');
        matrix_size = evalin('base','matrix_size');
        iFreq_raw = evalin('base','iFreq_raw');
        iMag = sqrt(sum(abs(iField).^2,4));
        assignin('base','iMag',iMag);
        set(handles.statusText,'String','Busy Performing Region Growth Unwrapping');
        drawnow
        iFreq = unwrapPhase(iMag,iFreq_raw,matrix_size);
        set(handles.statusText,'String','Done Performing Region Growth Unwrapping');
        assignin('base','iFreq',iFreq);
    otherwise 
        disp('fail');
end
catch
    FitButton_Callback(hObject, eventdata, handles);
    UnwrapButton_Callback(hObject, eventdata, handles);
end

% --- Executes on button press in RDFRB.
function RDFButton_Callback(hObject, eventdata, handles)
% hObject    handle to RDFRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(get(handles.BFRSel,'SelectedObject'),'Tag')
    case 'RDFRB'
        try
        iFreq = evalin('base','iFreq');
        N_std = evalin('base','N_std');
        Mask = evalin('base','Mask');
        matrix_size = evalin('base','matrix_size');
        voxel_size = evalin('base','voxel_size');
        B0_dir = evalin('base','B0_dir');
        set(handles.statusText,'String','Busy Performing BFR...');
        drawnow
        RDF = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);
        set(handles.statusText,'String','Done Performing BFR');
        assignin('base','RDF',RDF);
        catch ME2
            disp(['ME2' ME2.identifier]);
            UnwrapButton_Callback(hObject, eventdata, handles);
            RDFButton_Callback(hObject, eventdata, handles)
        end
    case 'LBVButton'
        try
        Mask = evalin('base','Mask');
        matrix_size = evalin('base','matrix_size');
        voxel_size = evalin('base','voxel_size');
        iFreq = evalin('base','iFreq');
        set(handles.statusText,'String','Busy Performing BFR...');
        drawnow
        RDF = LBV(iFreq,Mask,matrix_size,voxel_size);
        set(handles.statusText,'String','Done Performing BFR');
        assignin('base','RDF',RDF);
        catch ME3
            disp(['ME3' ME3.identifier]);
            UnwrapButton_Callback(hObject, eventdata, handles);
            RDFButton_Callback(hObject, eventdata, handles);
        end
    otherwise
        disp('Error');
end


% --- Executes on button press in MEDIButton.
function MEDIButton_Callback(hObject, eventdata, handles)
% hObject    handle to MEDIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(all(ismember(get(handles.SMVEdit,'String'), '0123456789+-.eEdD')))
    if(isempty(handles.varin))
        if(handles.SMVEnable)
            handles.varin = ['smv' ',' get(handles.SMVEdit,'String')];
            disp('SMV Entered');
        else
            disp('SMV Disabled');
        end
    end
%     disp(handles.varin);
end
if(all(ismember(get(handles.lambdaEditText,'String'), '0123456789+-.eEdD')))
    if(isempty(handles.varin))
        handles.varin = ['lambda' ',' get(handles.lambdaEditText,'String')];
    else
        handles.varin = ['lambda' ',' get(handles.lambdaEditText,'String') ',' handles.varin];
    end
%     disp(handles.varin);
end
if(all(ismember(get(handles.edgeEdit,'String'), '0123456789+-.eEdD')))
    if(isempty(handles.varin))
        handles.varin = ['percentage' ',' get(handles.edgeEdit,'String')];
    else
        handles.varin = ['percentage' ',' get(handles.edgeEdit,'String') ',' handles.varin];
    end
%     disp(handles.varin);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.varin = [handles.varin ',''filename'',' handles.filedir];
% disp(handles.varin);
try
    set(handles.statusText,'String','Busy Performing MEDI...');
    drawnow
    QSM = MEDI_L1(handles.varin);
    set(handles.statusText,'String','Done Performing MEDI');
    assignin('base','QSM',QSM);
catch ME
     disp(ME.identifier);
     saveRDFButton_Callback(hObject, eventdata, handles);
     MEDIButton_Callback(hObject, eventdata, handles);
end
    


% --- Executes on button press in Visu3DButton.
function Visu3DButton_Callback(hObject, eventdata, handles)
% hObject    handle to Visu3DButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
QSM = evalin('base','QSM');
voxel_size = evalin('base','voxel_size');
Mask = evalin('base','Mask');
Visu3D( QSM.*Mask, 'dimension',voxel_size);


% --- Executes when selected object is changed in uipanel9.
function uipanel9_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel9 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'LaplacianRadioButton'
        handles.unwrapSelection = get(eventdata.NewValue,'Tag');
    case 'RGRadioButton'
        handles.unwrapSelection = get(eventdata.NewValue,'Tag');
    otherwise
        disp(':(');
end


% --- Executes during object creation, after setting all properties.
function uipanel9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function lambdaEditText_Callback(hObject, eventdata, handles)
% hObject    handle to lambdaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambdaEditText as text
%        str2double(get(hObject,'String')) returns contents of lambdaEditText as a double


% --- Executes during object creation, after setting all properties.
function lambdaEditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lambdaEditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edgeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to edgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edgeEdit as text
%        str2double(get(hObject,'String')) returns contents of edgeEdit as a double


% --- Executes during object creation, after setting all properties.
function edgeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edgeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SMVEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SMVEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SMVEdit as text
%        str2double(get(hObject,'String')) returns contents of SMVEdit as a double


% --- Executes during object creation, after setting all properties.
function SMVEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SMVEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes when selected object is changed in uipanel13.
function uipanel13_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel13 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'defRB'
        set(handles.FitButton,'Enable','off');
        set(handles.UnwrapButton,'Enable','off');
        set(handles.RDFButton,'Enable','off');
        set(handles.MEDIButton,'Enable','off');
        set(handles.LaplacianRadioButton,'Enable','off');
        set(handles.RGRadioButton,'Enable','off');
        set(handles.RDFRB,'Enable','off');
        set(handles.LBVButton,'Enable','off');
        set(handles.SMVEdit,'Enable','off');
        set(handles.lambdaEditText,'Enable','off');
        set(handles.edgeEdit,'Enable','off');
        set(handles.saveRDFButton,'Enable','off');
        set(handles.buttonJust,'Enable','on');
    case 'advRB'
        set(handles.FitButton,'Enable','on');
        set(handles.UnwrapButton,'Enable','on');
        set(handles.RDFButton,'Enable','on');
        set(handles.MEDIButton,'Enable','on');
        set(handles.LaplacianRadioButton,'Enable','on');
        set(handles.RGRadioButton,'Enable','on');
        set(handles.RDFRB,'Enable','on');
        set(handles.LBVButton,'Enable','on');
        set(handles.SMVEdit,'Enable','on');
        set(handles.lambdaEditText,'Enable','on');
        set(handles.edgeEdit,'Enable','on');
        set(handles.saveRDFButton,'Enable','on');
        set(handles.buttonJust,'Enable','off');
    otherwise
        disp(':(');
end


% --- Executes on button press in buttonJust.
function buttonJust_Callback(hObject, eventdata, handles)
% hObject    handle to buttonJust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.statusText,'String','Busy Loading Data...');
    drawnow
[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir] = Read_DICOM(get(handles.DataLocationEditable,'string'));
set(handles.statusText,'String','Busy Loading Data...');
% xxx=GE or Siemens or Philips. 
% Other formats may be supported in the future
iMag = sqrt(sum(abs(iField).^2,4));
set(handles.statusText,'String','Busy Fitting ppm...');
    drawnow
[iFreq_raw N_std] = Fit_ppm_complex(iField);
set(handles.statusText,'String','Done Fitting ppm');
%%%%%Estimate the frequency offset in each of the voxel using a 
%%%%%complex fitting %%%%
set(handles.statusText,'String','Busy Unwrapping...');
drawnow
iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
set(handles.statusText,'String','Done Unwrapping');
% Spatial phase unwrapping %%%
% if large fringe lines persists, try 
 %iFreq = unwrapLaplacian(iFreq_raw, matrix_size, voxel_size);
% Mask = genMask(iField, voxel_size);
Mask = BET(iMag, matrix_size, voxel_size);
R2s = arlo(TE, abs(iField));
Mask_CSF = extract_CSF(R2s, Mask, voxel_size);
set(handles.statusText,'String','Busy Performing BF Removal...');
drawnow
RDF = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);
set(handles.statusText,'String','Done BF Removal');
%%%% Background field removal using Projection onto Dipole Fields
%%%% NMR Biomed 2011;24(9):1129-36.
%%%% MRM 2010;63(1):194-206
% Before running MEDI, variables need to be saved
save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
     voxel_size delta_TE CF B0_dir Mask_CSF;
set(handles.statusText,'String','Busy Performing MEDI...');
drawnow
QSM = MEDI_L1('lambda',1000,'percentage',0.9, 'smv',5);
set(handles.statusText,'String','Done MEDI');
assignin('base','QSM',QSM);
% morphology enabled dipole inversion


% --- Executes on button press in LoadMaskButton.
function LoadMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FilePath = get(handles.MaskLocationEditable,'string');
[pathstr,name,ext] = fileparts(FilePath);

switch ext
    case '.img'
        disp('Loading Mask');
        matrix_size = evalin('base','matrix_size');
        fid = fopen(FilePath);
        Mask = fread(fid,inf,'ushort');
        fclose(fid);
        Mask = reshape(Mask,matrix_size);
        assignin('base','Mask',Mask);
        disp('Done');
    case '.mat'
        disp('Loading Mask');
        mmask = load(FilePath,'-mat');
        Mask = mmask.Mask;
        assignin('base','Mask',Mask);
        disp('Done');
    otherwise 
        disp('invalid');
end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in updateLBButton.
function updateLBButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateLBButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
vars = evalin('base','who');
set(handles.listbox1,'String',vars);


% --- Executes on button press in visButton.
function visButton_Callback(hObject, eventdata, handles)
% hObject    handle to visButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varName = get(handles.listbox1,'String');
varValue = get(handles.listbox1,'Value');
varrrr = varName(varValue);
varrr = cell2mat(varrrr);
afk = evalin('base',varrr);
vis(afk);


% --- Executes on button press in saveRDFButton.
function saveRDFButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveRDFButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
handles.resultFolder = [get(handles.DataLocationEditable,'string') '_result'];
mkdir(handles.resultFolder);
iFreq = evalin('base','iFreq');
N_std = evalin('base','N_std');
Mask = evalin('base','Mask');
try
    Mask_CSF = evalin('base','Mask_CSF');
catch
    Mask_CSF = [];
end
matrix_size = evalin('base','matrix_size');
voxel_size = evalin('base','voxel_size');
B0_dir = evalin('base','B0_dir');
RDF = evalin('base','RDF');
iFreq_raw = evalin('base','iFreq_raw');
iMag = evalin('base','iMag');
delta_TE = evalin('base','delta_TE');
CF = evalin('base','CF');
handles.filedir = fullfile(handles.resultFolder, 'RDF.mat');
save (handles.filedir, 'RDF', 'iFreq', 'iFreq_raw', 'iMag', 'N_std', 'Mask', 'matrix_size',...
     'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'Mask_CSF');
 addpath(handles.resultFolder);
 guidata(hObject, handles);
catch ME1
    disp(ME1.identifier);
    RDFButton_Callback(hObject, eventdata, handles);

end


% --- Executes on button press in saveDICOMButton.
function saveDICOMButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveDICOMButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saveDir = uigetdir();
if (saveDir == 0)
    return
end
QSM = evalin('base','QSM');
write_QSM_dir(QSM,get(handles.DataLocationEditable,'string'),saveDir);


% --- Executes on button press in checkbox1.

function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
handles.SMVEnable = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function statusText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to statusText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in StripButton.
function StripButton_Callback(hObject, eventdata, handles)
% hObject    handle to StripButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mask = evalin('base','Mask');
Mask = stripBD(Mask,1);
assignin('base','Mask',Mask);
set(handles.statusText,'String','Mask Stripped');



function MultiEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MultiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MultiEdit as text
%        str2double(get(hObject,'String')) returns contents of MultiEdit as a double


% --- Executes during object creation, after setting all properties.
function MultiEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MultiEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MultiBrowse.
function MultiBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to MultiBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
paths= uigetdir();

if(isempty(get(handles.listbox2,'String')))
    new_name = [{paths}];
    set(handles.listbox2,'String',new_name);
    assignin('base','new_name',new_name);
else
    initial_name = cellstr(get(handles.listbox2,'String'));
    new_name = [initial_name;{paths}];
    set(handles.listbox2,'String',new_name);
    assignin('base','new_name',new_name);
end


% --- Executes on button press in MultiProcess.
function MultiProcess_Callback(hObject, eventdata, handles)
% hObject    handle to MultiProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_name = evalin('base','new_name');
[x,y]=size(new_name);
for i = 1:x
    X = sprintf('%d of %d is being processed.',i,x);
    disp(X);
    FilePath = new_name{i,1};
    set(handles.DataLocationEditable,'string',FilePath);
    set(handles.statusText,'String','Busy Loading Data... /n Please Wait');
    drawnow
    [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir] = Read_DICOM(FilePath);
    set(handles.statusText,'String','Done Loading Data');
    handles.loadedDataFolder = FilePath;
    assignin('base','iField',iField);
    assignin('base','voxel_size',voxel_size);
    assignin('base','matrix_size',matrix_size);
    assignin('base','CF',CF);
    assignin('base','delta_TE',delta_TE);
    assignin('base','TE',TE);
    assignin('base','B0_dir',B0_dir);
    iField = evalin('base','iField');
    voxel_size = evalin('base','voxel_size');
%     Mask = genMask(iField, voxel_size);
    iMag = sqrt(sum(abs(iField).^2,4));
    assignin('base','iMag',iMag);
    Mask = BET(iMag, matrix_size, voxel_size);
    R2s = arlo(TE, abs(iField));
    Mask_CSF = extract_CSF(R2s, Mask, voxel_size);
    assignin('base','Mask',Mask);
    assignin('base','R2s',R2s);
    assignin('base','Mask_CSF',Mask_CSF);
    set(handles.MaskLocationEditable,'Enable','off');
    set(handles.BrowseMaskButton,'Enable','off');
    set(handles.LoadMaskButton,'Enable','off');
    set(handles.statusText,'String','Mask Automatically Generated');
    drawnow
    FitButton_Callback(hObject, eventdata, handles);
    UnwrapButton_Callback(hObject, eventdata, handles);
    RDFButton_Callback(hObject, eventdata, handles);
    saveRDFButton_Callback(hObject, eventdata, handles);
    MEDIButton_Callback(hObject, eventdata, handles);
    QSM = evalin('base','QSM');
    write_QSM_dir(QSM,get(handles.DataLocationEditable,'string'),handles.resultFolder);
    
    
    guidata(hObject, handles);
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearDataDirButton.
function clearDataDirButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearDataDirButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.listbox2, 'String', '');
evalin('base', 'clearvars new_name');


% --- Executes during object creation, after setting all properties.
function uipanel13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
