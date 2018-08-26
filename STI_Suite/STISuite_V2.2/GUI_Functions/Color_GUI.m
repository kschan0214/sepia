function varargout = Color_GUI(varargin)
% Wei Li, Duke University, May 2010
% Wei Li, Updated on 10/28/2010 to include 4D display functionality
%


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Color_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Color_GUI_OutputFcn, ...
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


% --- Executes just before Color_GUI is made visible.
function Color_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Color_GUI (see VARARGIN)

% Choose default command line output for Color_GUI

handles.output = hObject;
handles.currentContour=1;
try
    set(handles.Load3DRgbEdit,'string',varargin{1});
    Visualize(hObject,handles,1)
end
% UIWAIT makes Color_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Color_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function LoadImages_Callback(hObject, eventdata, handles)
LoadImagesRGB(hObject,handles,1);

function Transforms_Callback(hObject, eventdata, handles)
Transforms(hObject,handles,1);


function LoadBinFIle_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,2);

function FrameSelect_Callback(hObject, eventdata, handles)

value = round(get(hObject,'Value'));
set(hObject,'Value',value);
set(handles.Frame,'String',num2str(value));
SliceSelection(hObject,handles);

function FrameSelect_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Frame_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FrameSelect3D_ButtonDownFcn(hObject, eventdata, handles)
SliceSelection(hObject,handles);

function SaveMag2Nii_Callback(hObject, eventdata, handles)
SaveMag2Nii(hObject,handles,1);

function Intensity_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,1);

function Intensity_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function Contrast_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,2);

function Contrast_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Load3DRgbEdit_Callback(hObject, eventdata, handles)
%set(handles.nDims,'String','Intensity');

VisualizeRGB(hObject,handles,1);

function VariableSelect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AutoScale1_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,3)

function Save2MatFile_Callback(hObject, eventdata, handles)
Save2MatFile(hObject,handles);

function ZeroPad_Callback(hObject, eventdata, handles)
ZeroPad(hObject,handles,1);

function CalculatePhase_Callback(hObject, eventdata, handles)
Transforms(hObject,handles,3);

function RestoreOrigionalSize_Callback(hObject, eventdata, handles)
ZeroPad(hObject,handles,2);

function Assign2Base_Callback(hObject, eventdata, handles)
Assign2Base(hObject,handles,1);

function Assign2Base_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function CalcMagn_Callback(hObject, eventdata, handles)
Transforms(hObject,handles,2);

function LoadNiiFile_Callback(hObject, eventdata, handles)
Normalization(hObject,handles,3);

function Load2ndImag_Callback(hObject, eventdata, handles)
LoadImagesRGB(hObject,handles,2);




function RightImageLoad_Callback(hObject, eventdata, handles)
VisualizeRGB(hObject,handles,2);

function RightImageLoad_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RightImageSave_Callback(hObject, eventdata, handles)
Assign2Base(hObject,handles,2);

function RightImageSave_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FrameSelect3D_Callback(hObject, eventdata, handles)

value = floor(get(hObject,'Value'));
set(handles.Frame,'String',num2str(value));

SliceSelectionColorGUI(hObject,handles);

function FrameSelect3D_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function Load3DRgbEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in FlipUpDown.
function FlipUpDown_Callback(hObject, eventdata, handles)
SliceSelectionColorGUI(hObject,handles);



% --- Executes on button press in Rotate90.
function Rotate90_Callback(hObject, eventdata, handles)
SliceSelectionColorGUI(hObject,handles);

% --- Executes on slider movement.
function Coilslider_Callback(hObject, eventdata, handles)

value = round(get(hObject,'Value'));
set(hObject,'Value',value);
set(handles.CoilNum,'String',num2str(value));
SliceSelection(hObject,handles);



% --- Executes during object creation, after setting all properties.
function Coilslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on FrameSelect3D and none of its controls.
function FrameSelect3D_KeyPressFcn(hObject, eventdata, handles)

switch eventdata.Character
    case '4'
        value = floor(get(handles.FrameSelect3D,'Value'));
        if  value <2
        else
            set(handles.FrameSelect3D,'Value',value-1);
            set(handles.Frame,'String',num2str(value-1));
            SliceSelection(hObject,handles);
        end
        
    case '6'
        value = floor(get(handles.FrameSelect3D,'Value'));
        if  value > get(handles.FrameSelect3D,'max')-1
        else
            set(handles.FrameSelect3D,'Value',value+1);
            set(handles.Frame,'String',num2str(value+1));
           SliceSelection(hObject,handles);
        end
    case '8' 
        value = floor(get(handles.Coilslider,'Value'));
        if  value > get(handles.Coilslider,'max')-1
        else
            set(handles.Coilslider,'Value',value+1);
            set(handles.CoilNum,'String',num2str(value+1));
            SliceSelection(hObject,handles);
        end
    case '2'
        value = floor(get(handles.Coilslider,'Value'));
        if value<2 
        else
            set(handles.Coilslider,'Value',value-1);
            set(handles.CoilNum,'String',num2str(value-1));
            SliceSelection(hObject,handles);
        end
        
end



function Load3DRgbEditRight_Callback(hObject, eventdata, handles)
VisualizeRGB(hObject,handles,2);

% --- Executes during object creation, after setting all properties.
function Load3DRgbEditRight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function RGBlow_Callback(hObject, eventdata, handles)
SliceSelectionColorGUI(hObject,handles);



% --- Executes during object creation, after setting all properties.
function RGBlow_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RGBhigh_Callback(hObject, eventdata, handles)
SliceSelectionColorGUI(hObject,handles);


% --- Executes during object creation, after setting all properties.
function RGBhigh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Frame_Callback(hObject, eventdata, handles)
framenumber=str2double(get(handles.Frame,'string'));
set(handles.FrameSelect3D,'value',framenumber);
SliceSelection(hObject,handles);


function RightLow_Callback(hObject, eventdata, handles)
SliceSelectionColorGUI(hObject,handles);


% --- Executes during object creation, after setting all properties.
function RightLow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RightHigh_Callback(hObject, eventdata, handles)
SliceSelectionColorGUI(hObject,handles);

% --- Executes during object creation, after setting all properties.
function RightHigh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press with focus on Coilslider and none of its controls.
function Coilslider_KeyPressFcn(hObject, eventdata, handles)


% --- Executes when selected object is changed in ImagePlane.
function ImagePlane_SelectionChangeFcn(hObject, eventdata, handles)

switch get(hObject,'string')
    case 'Sagittal (Y-Z)'
        handles.ImageAxisNum=1;
    case 'Coronal (X-Z)'
        handles.ImageAxisNum=2;
    case 'Axial (X-Y)'
        handles.ImageAxisNum=3;
end
guidata(hObject, handles);       

SliceSelectionColorGUI(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ImagePlane_CreateFcn(hObject, eventdata, handles)
handles.ImageAxisNum=3;
guidata(hObject, handles);       


% --- Executes on button press in SaveMag2Nii2.
function SaveMag2Nii2_Callback(hObject, eventdata, handles)
SaveMag2Nii(hObject,handles,2);
% hObject    handle to SaveMag2Nii2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LoadBinFIle2.
function LoadBinFIle2_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,4);

% hObject    handle to LoadBinFIle2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in readmrdata.
function readmrdata_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,5);



% --- Executes on button press in loadphase.
function loadphase_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,6);



% --- Executes when selected object is changed in uipanel10.
function uipanel10_SelectionChangeFcn(hObject, eventdata, handles)

set(handles.Coilslider, 'value', 1); 
switch get(hObject,'string')
    case 'PhiFiltered'
        handles.M2=handles.PhiFiltered;
    case 'UpdatedMask'
        handles.M2=handles.FinalMask;
    case 'X'
        handles.M2=handles.X;
end
guidata(hObject, handles);       

if strcmp(get(handles.nDims,'string'),'Intensity')
    SliceSelection(hObject,handles);
else
    SliceSelectionColorGUI(hObject,handles);
end
AdjustIntensity(hObject,handles,3)




% --------------------------------------------------------------------
function uitoggletool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Scaling.
% function Scaling_Callback(hObject, eventdata, handles)

% SliceSelection(hObject,handles);


% --------------------------------------------------------------------
function uitoggletool3_OffCallback(hObject, eventdata, handles)
set(handles.ImageAxes2,'XLim',get(handles.ImageAxes,'XLim'),...
                       'YLim',get(handles.ImageAxes,'YLim'));


% --------------------------------------------------------------------
function uitoggletool5_OffCallback(hObject, eventdata, handles)
set(handles.ImageAxes2,'XLim',get(handles.ImageAxes,'XLim'),...
                       'YLim',get(handles.ImageAxes,'YLim'));




% --- Executes on button press in Scaling.
function Scaling_Callback(hObject, eventdata, handles)
SliceSelectionColorGUI(hObject,handles);



function LeftResolution_Callback(hObject, eventdata, handles)
% hObject    handle to LeftResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LeftResolution as text
%        str2double(get(hObject,'String')) returns contents of LeftResolution as a double


% --- Executes during object creation, after setting all properties.
function LeftResolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LeftResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RightResolution_Callback(hObject, eventdata, handles)
% hObject    handle to RightResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RightResolution as text
%        str2double(get(hObject,'String')) returns contents of RightResolution as a double


% --- Executes during object creation, after setting all properties.
function RightResolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RightResolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to RightLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RightLow as text
%        str2double(get(hObject,'String')) returns contents of RightLow as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RightLow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to RightHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RightHigh as text
%        str2double(get(hObject,'String')) returns contents of RightHigh as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RightHigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
