function varargout = Gray_ROI(varargin)
% Wei Li, Duke University, May 2010
% Wei Li, Updated on 10/28/2010 to include 4D display functionality
%


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gray_ROI_OpeningFcn, ...
                   'gui_OutputFcn',  @Gray_ROI_OutputFcn, ...
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


% --- Executes just before Gray_ROI is made visible.
function Gray_ROI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gray_ROI (see VARARGIN)

% Choose default command line output for Gray_ROI

handles.output = hObject;
handles.currentContour=1;
guidata(hObject, handles);


% UIWAIT makes Gray_ROI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Gray_ROI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function LoadImages_Callback(hObject, eventdata, handles)
Brain_ROIs(hObject,handles);

function Transforms_Callback(hObject, eventdata, handles)
Transforms(hObject,handles,1);

function SavePhase2Bin_Callback(hObject, eventdata, handles)
SavePhase2Bin(hObject,handles)

function LoadBinFIle_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,2);

function FrameSelect_Callback(hObject, eventdata, handles)

value = floor(get(hObject,'Value'));
set(hObject,'Value',value);


set(handles.Frame,'String',num2str(value));
SliceSelectionROI(hObject,handles);

function FrameSelect_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Frame_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FrameSelect_ButtonDownFcn(hObject, eventdata, handles)
SliceSelectionROI(hObject,handles);

function SaveMag2Nii_Callback(hObject, eventdata, handles)
SaveMag2Nii(hObject,handles);

function Intensity_Callback(hObject, eventdata, handles)
AdjuGray_ROIntensity(hObject,handles,1);

function Intensity_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Contrast_Callback(hObject, eventdata, handles)
AdjuGray_ROIntensity(hObject,handles,2);

function Contrast_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function VariableSelect_Callback(hObject, eventdata, handles)
%set(handles.nDims,'String','Intensity');
VisualizeROI(hObject,handles,1);

function VariableSelect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Disp_SelectionChangeFcn(hObject, eventdata, handles)
SliceSelectionROI(hObject,handles);

function ColorScheme_SelectionChangeFcn(hObject, eventdata, handles)
ColorScheme(hObject,handles);

function AutoScale_Callback(hObject, eventdata, handles)
AdjustIntensityROI(hObject,handles,3)

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
LoadBinFile(hObject,handles,3);

function myZoom_Callback(hObject, eventdata, handles)


function RightImageLoad_Callback(hObject, eventdata, handles)
Assign2BaseROI(hObject,handles,1);
% value = round(get(hObject,'Value'));
% set(hObject,'Value',value);
% set(handles.CoilNum,'String',num2str(value));
% SliceSelectionROI(hObject,handles);


function RightImageLoad_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RightImageSave_Callback(hObject, eventdata, handles)
set(handles.status,'string',[])


function RightImageSave_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FrameSelect3D_Callback(hObject, eventdata, handles)

value = floor(get(hObject,'Value'));
set(handles.Frame,'String',num2str(value));

SliceSelection3D(hObject,handles);

function FrameSelect3D_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function Load3DRgbEdit_Callback(hObject, eventdata, handles)
set(handles.nDims,'String','RGB');

VisualizeRGB(hObject,handles,1);

function Load3DRgbEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AdjustContrastButton_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,5);


% --- Executes on button press in FollowLightContrast.
function FollowLightContrast_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,6);



% --- Executes on button press in AdjustContrastR.
function AdjustContrastR_Callback(hObject, eventdata, handles)
AdjustIntensityROI(hObject,handles,7);


% --- Executes on button press in FollowRightContrast.
function FollowRightContrast_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,8);

% --- Executes when selected object is changed in FlipChange.

% --- Executes on button press in FlipUpDown.
function FlipUpDown_Callback(hObject, eventdata, handles)
SliceSelectionROI(hObject,handles);


% --- Executes on button press in Rotate90.
function Rotate90_Callback(hObject, eventdata, handles)
SliceSelectionROI(hObject,handles);

% --- Executes on slider movement.
function Coilslider_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));

set(hObject,'Value',value);
set(handles.CoilNum,'String',num2str(value));
%guidata(gcbo,handles);
%ROI_Update(hObject,handles);
try
SliceSelectionROI(hObject,handles);
end
AdjustIntensityROI(hObject,handles,3)
% --- Executes during object creation, after setting all properties.
function Coilslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on FrameSelect and none of its controls.
function FrameSelect_KeyPressFcn(hObject, eventdata, handles)

switch eventdata.Character
    case '4'
        value = floor(get(handles.FrameSelect,'Value'));
        if  value <2
        else
            set(handles.FrameSelect,'Value',value-1);
            set(handles.Frame,'String',num2str(value-1));
            SliceSelectionROI(hObject,handles);
        end
        
    case '6'
        value = floor(get(handles.FrameSelect,'Value'));
        if  value > get(handles.FrameSelect,'max')-1
        else
            set(handles.FrameSelect,'Value',value+1);
            set(handles.Frame,'String',num2str(value+1));
           SliceSelectionROI(hObject,handles);
        end
    case '8' 
        value = floor(get(handles.Coilslider,'Value'));
        if  value > get(handles.Coilslider,'max')-1
        else
            set(handles.Coilslider,'Value',value+1);
            set(handles.CoilNum,'String',num2str(value+1));
            SliceSelectionROI(hObject,handles);
        end
    case '2'
        value = floor(get(handles.Coilslider,'Value'));
        if value<2 
        else
            set(handles.Coilslider,'Value',value-1);
            set(handles.CoilNum,'String',num2str(value-1));
            SliceSelectionROI(hObject,handles);
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
SliceSelection3D(hObject,handles);



% --- Executes during object creation, after setting all properties.
function RGBlow_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ROIValue_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ROIValue_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Frame_Callback(hObject, eventdata, handles)
framenumber=str2double(get(handles.Frame,'string'));
set(handles.FrameSelect,'value',framenumber);
SliceSelectionROI(hObject,handles);


function RightLow_Callback(hObject, eventdata, handles)
SliceSelection3D(hObject,handles);


% --- Executes during object creation, after setting all properties.
function RightLow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RightHigh_Callback(hObject, eventdata, handles)
SliceSelection3D(hObject,handles);

% --- Executes during object creation, after setting all properties.
function RightHigh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on key press with focus on Coilslider and none of its controls.
function Coilslider_KeyPressFcn(hObject, eventdata, handles)



% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
value = floor(get(handles.FrameSelect,'value'));
if eventdata.VerticalScrollCount==1
    value=value+1;
else
    value=value-1;
end
if value <1
    value =1;
end
SS=size(handles.Mi);


if value >SS(handles.ImageAxisNum)
    value =SS(handles.ImageAxisNum);
end

set(handles.FrameSelect,'value',value);
guidata(hObject, handles);
SliceSelectionROI(hObject,handles);



% --- Executes when selected object is changed in Numberselect.
function Numberselect_SelectionChangeFcn(hObject, eventdata, handles)
handles.NumberSelect=str2double(get(hObject,'string'));
guidata(hObject, handles);
Brain_ROIs(hObject,handles)



% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
value = floor(get(handles.FrameSelect,'value'));
%eventdata.Key
switch eventdata.Key
    case 'uparrow'
        value=value+1;
set(handles.FrameSelect,'value',value);
    case 'downarrow'
        value=value-1;
set(handles.FrameSelect,'value',value);
    case 'w'
        value=value+1;
 set(handles.FrameSelect,'value',value);
   case 's'
        value=value-1;
    set(handles.FrameSelect,'value',value);
    case 'hyphen'
    value=(get(handles.ROISize,'Value'))-1;
    if value <1; value=1;end
    if value>5; value=5;end
    set(handles.ROISize,'Value',value);
    set(handles.radius0,'String',num2str(value));
    case 'equal'
    value=(get(handles.ROISize,'Value'))+1;
      if value <1; value=1;end
    if value>5; value=5;end
    set(handles.ROISize,'Value',value);
    set(handles.radius0,'String',num2str(value));

end
guidata(hObject, handles);
SliceSelectionROI(hObject,handles);


% --------------------------------------------------------------------
function uitoggletool3_OffCallback(hObject, eventdata, handles)
BrainZoom(hObject,handles)


% --------------------------------------------------------------------
function uitoggletool9_OffCallback(hObject, eventdata, handles)
BrainZoom(hObject,handles)



% --- Executes when selected object is changed in ImagePlane.
function ImagePlane_SelectionChangeFcn(hObject, eventdata, handles)

switch get(hObject,'string')
    case 'Axial (X-Y)'
        handles.ImageAxisNum=3;
    case 'Sagittal (Y-Z)'
        handles.ImageAxisNum=1;
    case 'Coronal (X-Z)'
        handles.ImageAxisNum=2;
end

guidata(hObject, handles);       

SliceSelectionROI(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ImagePlane_CreateFcn(hObject, eventdata, handles)

handles.ImageAxisNum=3;
guidata(hObject, handles);       




function radius0_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function radius0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveROI.
function SaveROI_Callback(hObject, eventdata, handles)
Assign2BaseROI_Gray(hObject,handles,2);



% --- Executes on slider movement.
function ROISize_Callback(hObject, eventdata, handles)
value = round(get(hObject,'Value'));
set(hObject,'Value',value);
set(handles.radius0,'String',num2str(value));



% --- Executes during object creation, after setting all properties.
function ROISize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROISize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
LoadBinFileROI(hObject,handles,5);

% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)

Assign2BaseROI_Gray(hObject,handles,3);

% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)

Assign2BaseROI_Gray(hObject,handles,4);

% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uitoggletool5_OffCallback(hObject, eventdata, handles)
BrainZoom(hObject,handles)
% hObject    handle to uitoggletool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
LoadBinFileROI(hObject,handles,1);
set(handles.Number1,'value',1);
set(handles.Number2,'value',0);
set(handles.Number3,'value',0);
set(handles.number4,'value',0);
set(handles.number5,'value',0);
set(handles.number6,'value',0);
set(handles.number7,'value',0);
set(handles.number8,'value',0);
set(handles.number9,'value',0);

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)

value = get(hObject,'Value');
a= size(get(handles.himage2,'AlphaData'));

if value==0
    set(handles.himage2,'AlphaData',ones(a));
else
    SliceSelectionROI(hObject,handles)
end




function edit17_Callback(hObject, eventdata, handles)
Assign3ROI(hObject,handles,1);


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function AlphaValueSlider_Callback(hObject, eventdata, handles)
SliceSelectionROI(hObject,handles);


% --- Executes during object creation, after setting all properties.
function AlphaValueSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AlphaValueSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function LeftResolution_Callback(hObject, eventdata, handles)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LeftResolution as text
%        str2double(get(hObject,'String')) returns contents of LeftResolution as a double


% --- Executes during object creation, after setting all properties.
function LeftResolution_CreateFcn(hObject, eventdata, handles)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Scaling.
function Scaling_Callback(hObject, eventdata, handles)
% hObject    handle to Scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Scaling


% --- Executes on button press in loadniiroi.
function loadniiroi_Callback(hObject, eventdata, handles)
Assign2BaseROI_Gray(hObject,handles,5);


% --- Executes on button press in saveniiroi.
function saveniiroi_Callback(hObject, eventdata, handles)
Assign2BaseROI_Gray(hObject,handles,6);



% --- Executes on button press in ClearLabels.
function ClearLabels_Callback(hObject, eventdata, handles)

handles.M2=ones(size(handles.Mi,1),size(handles.Mi,2),size(handles.Mi,3));
guidata(gcbo,handles);
try
SliceSelectionROI(hObject,handles);
end




function AddMoreImages_Callback(hObject, eventdata, handles)

VisualizeROI(hObject,handles,2);



% --- Executes during object creation, after setting all properties.
function AddMoreImages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AddMoreImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



