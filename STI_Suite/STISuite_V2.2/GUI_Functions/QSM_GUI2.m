function varargout = QSM_GUI2(varargin)
% Wei Li, Duke University, May 2010
% Wei Li, Updated on 10/28/2010 to include 4D display functionality
%


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @QSM_GUI2_OpeningFcn, ...
                   'gui_OutputFcn',  @QSM_GUI2_OutputFcn, ...
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


% --- Executes just before QSM_GUI2 is made visible.
function QSM_GUI2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to QSM_GUI2 (see VARARGIN)

% Choose default command line output for QSM_GUI2

handles.output = hObject;
handles.currentContour=1;
try
    set(handles.VariableSelect,'string',varargin{1});
    Visualize(hObject,handles,1)
end


% UIWAIT makes QSM_GUI2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
load IconImage
Image=IconImage;
handles.Mi=Image;
Image=permute(Image,[2 1 3]);
Image=flipdim(Image,1);
handles.himage1=imshow(Image(:,:,1),[],'parent',handles.ImageAxes);
set(handles.ImageAxes,'Ydir','reverse');

Image=BrainIcon;
handles.M2=Image;
Image=permute(Image,[2 1 3]);
Image=flipdim(Image,1);
handles.himage2=imshow(Image(:,:,1),[],'parent',handles.ImageAxes2);
set(handles.ImageAxes2,'Ydir','reverse');


AA=get(handles.QSMProcessing,'position');
x0=get(handles.PhaseProcessing,'position');
set(handles.QSMProcessing,'position',[x0(1) AA(2) AA(3) AA(4)]);

x0=get(handles.PhaseProcessing,'position');
AA=get(handles.VisualizeSelect,'position');
set(handles.VisualizeSelect,'position',[x0(1) AA(2) AA(3) AA(4)]);

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = QSM_GUI2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function LoadImages_Callback(hObject, eventdata, handles)
LoadImages(hObject,handles,1);

function Transforms_Callback(hObject, eventdata, handles)
Transforms(hObject,handles,1);

function Thres_Callback(hObject, eventdata, handles)
mythresh(hObject,handles);

function LoadBinFIle_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,22);

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

function FrameSelect_ButtonDownFcn(hObject, eventdata, handles)
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

function VariableSelect_Callback(hObject, eventdata, handles)
%set(handles.nDims,'String','Intensity');
Visualize(hObject,handles,1);

function VariableSelect_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Disp_SelectionChangeFcn(hObject, eventdata, handles)
SliceSelection(hObject,handles);

function ColorScheme_SelectionChangeFcn(hObject, eventdata, handles)
ColorScheme(hObject,handles);

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
LoadImages(hObject,handles,2);




function RightImageLoad_Callback(hObject, eventdata, handles)
Visualize(hObject,handles,2);

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
AdjustIntensity(hObject,handles,7);


% --- Executes on button press in FollowRightContrast.
function FollowRightContrast_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,8);

% --- Executes when selected object is changed in FlipChange.

% --- Executes on button press in FlipUpDown.
function FlipUpDown_Callback(hObject, eventdata, handles)
SliceSelection(hObject,handles);



% --- Executes on button press in Rotate90.
function Rotate90_Callback(hObject, eventdata, handles)
SliceSelection(hObject,handles);

% --- Executes on slider movement.
function Coilslider_Callback(hObject, eventdata, handles)

if length(size(handles.Mi))==3
    set(handles.Coilslider, 'value', 1); 
    set(handles.CoilNum, 'string', '1'); 
else
    sss=size(handles.Mi,4);
    set(handles.Coilslider,'Min',1,'Max',sss,'SliderStep',[1/(sss-1) 2/(sss-1)]); 
end

value = round(get(hObject,'Value'));
set(hObject,'Value',value);
set(handles.CoilNum,'String',num2str(value));
SliceSelection(hObject,handles);



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
            SliceSelection(hObject,handles);
        end
        
    case '6'
        value = floor(get(handles.FrameSelect,'Value'));
        if  value > get(handles.FrameSelect,'max')-1
        else
            set(handles.FrameSelect,'Value',value+1);
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
SliceSelection3D(hObject,handles);



% --- Executes during object creation, after setting all properties.
function RGBlow_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RGBhigh_Callback(hObject, eventdata, handles)
SliceSelection3D(hObject,handles);


% --- Executes during object creation, after setting all properties.
function RGBhigh_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Frame_Callback(hObject, eventdata, handles)
framenumber=str2double(get(handles.Frame,'string'));
set(handles.FrameSelect,'value',framenumber);
SliceSelection(hObject,handles);


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


% --- Executes when selected object is changed in ImagePlane.
function ImagePlane_SelectionChangeFcn(hObject, eventdata, handles)

SS=size(handles.Mi);
switch get(hObject,'string')
    case 'Sagittal (Y-Z)'
        handles.ImageAxisNum=1;
        set(handles.FrameSelect,'value',round(SS(1)/2));
        set(handles.Frame,'String',num2str(round(SS(1)/2)));
    case 'Coronal (X-Z)'
        handles.ImageAxisNum=2;
        set(handles.FrameSelect,'value',round(SS(2)/2));
        set(handles.Frame,'String',num2str(round(SS(2)/2)));
    case 'Axial (X-Y)'
        handles.ImageAxisNum=3;
        set(handles.FrameSelect,'value',round(SS(3)/2));
        set(handles.Frame,'String',num2str(round(SS(3)/2)));
end
guidata(hObject, handles);       
SliceSelection(hObject,handles);


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
LoadBinFile(hObject,handles,44);

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
    SliceSelection3D(hObject,handles);
end
AdjustIntensity(hObject,handles,3)




% --------------------------------------------------------------------
function uitoggletool3_ClickedCallback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uitoggletool3_OffCallback(hObject, eventdata, handles)
set(handles.ImageAxes2,'XLim',get(handles.ImageAxes,'XLim'),...
                       'YLim',get(handles.ImageAxes,'YLim'));


% --------------------------------------------------------------------
function uitoggletool5_OffCallback(hObject, eventdata, handles)
set(handles.ImageAxes2,'XLim',get(handles.ImageAxes,'XLim'),...
                       'YLim',get(handles.ImageAxes,'YLim'));




% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
SliceSelection(hObject,handles);



function LeftResolution_Callback(hObject, eventdata, handles)

phaseres=get(handles.LeftResolution,'string');
set(handles.VoxelSizeForQSM,'string',phaseres);
set(handles.Voxelsize5,'string',phaseres);
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




% --- Executes on button press in AutoScale.
function AutoScale_Callback(hObject, eventdata, handles)
AdjustIntensity(hObject,handles,3)






function BrainMaskName_Callback(hObject, eventdata, handles)
Visualize(hObject,handles,6);
BrainMaskName=get(handles.BrainMaskName,'string');
set(handles.NewMaskName,'string',BrainMaskName);



% --- Executes during object creation, after setting all properties.
function BrainMaskName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BrainMaskName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Loadmask.
function Loadmask_Callback(hObject, eventdata, handles)
LoadImages(hObject,handles,3);



% hObject    handle to Loadmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function TE1_Callback(hObject, eventdata, handles)
TE1=str2double(get(handles.TE1,'string'));
DeltaTE=str2double(get(handles.DeltaTE,'string'));
SS=size(handles.Mi);
if length(SS)==3
   SS(4)=1; 
end
TEs=TE1+DeltaTE*(0:SS(4)-1);
set(handles.TE_value,'string',num2str(sum(TEs)))

% hObject    handle to TE1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TE1 as text
%        str2double(get(hObject,'String')) returns contents of TE1 as a double


% --- Executes during object creation, after setting all properties.
function TE1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TE1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function H_Vector_Callback(hObject, eventdata, handles)
% hObject    handle to H_Vector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of H_Vector as text
%        str2double(get(hObject,'String')) returns contents of H_Vector as a double


% --- Executes during object creation, after setting all properties.
function H_Vector_CreateFcn(hObject, eventdata, handles)
% hObject    handle to H_Vector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
QSM_Functions(hObject,handles,'VSHARP')

% --- Executes on button press in HARPERELLA.
function HARPERELLA_Callback(hObject, eventdata, handles)
CalliHARPERELLA(hObject,handles,'HARPERELLA')

set(handles.VisualizeSelect,'visible','off')
set(handles.PhaseProcessing,'visible','off')
set(handles.QSMProcessing,'visible','on')
set(handles.QSMSelect,'value',1)
drawnow;

function ProcessedPhaseName_Callback(hObject, eventdata, handles)
Visualize(hObject,handles,3);


% --- Executes during object creation, after setting all properties.
function ProcessedPhaseName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProcessedPhaseName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NewMaskName_Callback(hObject, eventdata, handles)
Visualize(hObject,handles,4);


% --- Executes during object creation, after setting all properties.
function NewMaskName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewMaskName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PadSizeForQSM_Callback(hObject, eventdata, handles)
% hObject    handle to PadSizeForQSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PadSizeForQSM as text
%        str2double(get(hObject,'String')) returns contents of PadSizeForQSM as a double


% --- Executes during object creation, after setting all properties.
function PadSizeForQSM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PadSizeForQSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VoxelSizeForQSM_Callback(hObject, eventdata, handles)
% hObject    handle to VoxelSizeForQSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VoxelSizeForQSM as text
%        str2double(get(hObject,'String')) returns contents of VoxelSizeForQSM as a double


% --- Executes during object creation, after setting all properties.
function VoxelSizeForQSM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VoxelSizeForQSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in QSM_LSQR.
function QSM_LSQR_Callback(hObject, eventdata, handles)
QSM_LSQR_Callback_Fcn(hObject,handles,'QSM')


function B0Value_Callback(hObject, eventdata, handles)
% hObject    handle to B0Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of B0Value as text
%        str2double(get(hObject,'String')) returns contents of B0Value as a double


% --- Executes during object creation, after setting all properties.
function B0Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to B0Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TE_value_Callback(hObject, eventdata, handles)
% hObject    handle to TE_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TE_value as text
%        str2double(get(hObject,'String')) returns contents of TE_value as a double


% --- Executes during object creation, after setting all properties.
function TE_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TE_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function UniqueID_Callback(hObject, eventdata, handles)
% hObject    handle to UniqueID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UniqueID as text
%        str2double(get(hObject,'String')) returns contents of UniqueID as a double


% --- Executes during object creation, after setting all properties.
function UniqueID_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UniqueID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Niter_Callback(hObject, eventdata, handles)
% hObject    handle to Niter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Niter as text
%        str2double(get(hObject,'String')) returns contents of Niter as a double


% --- Executes during object creation, after setting all properties.
function Niter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Niter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PhaseIter_Callback(hObject, eventdata, handles)
% hObject    handle to PhaseIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PhaseIter as text
%        str2double(get(hObject,'String')) returns contents of PhaseIter as a double


% --- Executes during object creation, after setting all properties.
function PhaseIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhaseIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function phaseRadius_Callback(hObject, eventdata, handles)
% hObject    handle to phaseRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phaseRadius as text
%        str2double(get(hObject,'String')) returns contents of phaseRadius as a double


% --- Executes during object creation, after setting all properties.
function phaseRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function VSHARPRadius_Callback(hObject, eventdata, handles)
% hObject    handle to VSHARPRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VSHARPRadius as text
%        str2double(get(hObject,'String')) returns contents of VSHARPRadius as a double


% --- Executes during object creation, after setting all properties.
function VSHARPRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VSHARPRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectMaskNii.
function SelectMaskNii_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,45);



% --- Executes on button press in Deconvolution.
function Deconvolution_Callback(hObject, eventdata, handles)
% hObject    handle to Deconvolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Deconvolution


% --- Executes on button press in FastQSM.
function FastQSM_Callback(hObject, eventdata, handles)
FastQSMCallBack(hObject,handles,'EstimateX0')



% --- Executes on button press in QSMKIWI.
function QSMKIWI_Callback(hObject, eventdata, handles)
% hObject    handle to QSMKIWI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
QSM_Functions(hObject,handles,'QSM_KIWI')



function cropsize_Callback(hObject, eventdata, handles)
% hObject    handle to cropsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cropsize as text
%        str2double(get(hObject,'String')) returns contents of cropsize as a double


% --- Executes during object creation, after setting all properties.
function cropsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tol_step1_Callback(hObject, eventdata, handles)
% hObject    handle to tol_step1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tol_step1 as text
%        str2double(get(hObject,'String')) returns contents of tol_step1 as a double


% --- Executes during object creation, after setting all properties.
function tol_step1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tol_step1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tol_step2_Callback(hObject, eventdata, handles)
% hObject    handle to tol_step2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tol_step2 as text
%        str2double(get(hObject,'String')) returns contents of tol_step2 as a double


% --- Executes during object creation, after setting all properties.
function tol_step2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tol_step2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fraction_kspace_Callback(hObject, eventdata, handles)
% hObject    handle to fraction_kspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fraction_kspace as text
%        str2double(get(hObject,'String')) returns contents of fraction_kspace as a double


% --- Executes during object creation, after setting all properties.
function fraction_kspace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fraction_kspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function YZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in phaseuwrapping.
function phaseuwrapping_Callback(hObject, eventdata, handles)
phaseuwrapping_CallbackFcn(hObject,handles,'Unwrapping')



function unwrappedphase_Callback(hObject, eventdata, handles)
Visualize(hObject,handles,5);



% --- Executes during object creation, after setting all properties.
function unwrappedphase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to unwrappedphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Customize_Callback(hObject, eventdata, handles)
guide QSM_GUI
% hObject    handle to Customize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Customization_Callback(hObject, eventdata, handles)

% hObject    handle to Customization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RunVSharp.
function RunVSharp_Callback(hObject, eventdata, handles)

VSHARPCallBackFcn(hObject,handles,'VSHARP')
set(handles.PhaseProcessing,'visible','off')
set(handles.QSMProcessing,'visible','on')
set(handles.VisualizeSelect,'visible','off')
set(handles.QSMSelect,'value',1)
drawnow;


% --- Executes on button press in pushbutton41.
function pushbutton41_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton41 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit43_Callback(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit43 as text
%        str2double(get(hObject,'String')) returns contents of edit43 as a double


% --- Executes during object creation, after setting all properties.
function edit43_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit44 as text
%        str2double(get(hObject,'String')) returns contents of edit44 as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton42.
function pushbutton42_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton42 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit45_Callback(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit45 as text
%        str2double(get(hObject,'String')) returns contents of edit45 as a double


% --- Executes during object creation, after setting all properties.
function edit45_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton43.
function pushbutton43_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton43 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton44.
function pushbutton44_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton44 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton45.
function pushbutton45_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton45 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton46.
function pushbutton46_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit46_Callback(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit46 as text
%        str2double(get(hObject,'String')) returns contents of edit46 as a double


% --- Executes during object creation, after setting all properties.
function edit46_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit46 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit47_Callback(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit47 as text
%        str2double(get(hObject,'String')) returns contents of edit47 as a double


% --- Executes during object creation, after setting all properties.
function edit47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton47.
function pushbutton47_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit48_Callback(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit48 as text
%        str2double(get(hObject,'String')) returns contents of edit48 as a double


% --- Executes during object creation, after setting all properties.
function edit48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit49_Callback(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit49 as text
%        str2double(get(hObject,'String')) returns contents of edit49 as a double


% --- Executes during object creation, after setting all properties.
function edit49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit50_Callback(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit50 as text
%        str2double(get(hObject,'String')) returns contents of edit50 as a double


% --- Executes during object creation, after setting all properties.
function edit50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function LoadRealData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadRealData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LoadRealData as text
%        str2double(get(hObject,'String')) returns contents of LoadRealData as a double


% --- Executes during object creation, after setting all properties.
function LoadRealData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadRealData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LoadImagData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImagData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LoadImagData as text
%        str2double(get(hObject,'String')) returns contents of LoadImagData as a double


% --- Executes during object creation, after setting all properties.
function LoadImagData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadImagData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadimagnii.
function loadimagnii_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,26);



% --- Executes on button press in loadimag.
function loadimag_Callback(hObject, eventdata, handles)
LoadImages(hObject,handles,7);


% --- Executes on button press in loadreal.
function loadreal_Callback(hObject, eventdata, handles)
LoadImages(hObject,handles,6);



% --- Executes on button press in loadrealnii.
function loadrealnii_Callback(hObject, eventdata, handles)
LoadBinFile(hObject,handles,23);




% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function DeltaTE_Callback(hObject, eventdata, handles)
TE1=str2double(get(handles.TE1,'string'));
DeltaTE=str2double(get(handles.DeltaTE,'string'));
SS=size(handles.Mi);
if length(SS)==3
   SS(4)=1; 
end
TEs=TE1+DeltaTE*(0:SS(4)-1);
set(handles.TE_value,'string',num2str(sum(TEs)));


% --- Executes during object creation, after setting all properties.
function DeltaTE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DeltaTE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in PanalSelection.
function PanalSelection_SelectionChangeFcn(hObject, eventdata, handles)
set(handles.CoilNum,'String',1);
switch get (hObject,'string')
    case 'phase'
        set(handles.PhaseProcessing,'visible','on')
        set(handles.QSMProcessing,'visible','off')
        set(handles.VisualizeSelect,'visible','off')
    case 'QSM'
        set(handles.PhaseProcessing,'visible','off')
        set(handles.QSMProcessing,'visible','on')
        set(handles.VisualizeSelect,'visible','off')
        
    case 'Visualize'
        set(handles.PhaseProcessing,'visible','off')
        set(handles.QSMProcessing,'visible','off')
        set(handles.VisualizeSelect,'visible','on')
end

function VisualizeSelect_SelectionChangeFcn(hObject, eventdata, handles)


function CombinedLaplacian_Callback(hObject, eventdata, handles)
CalcCombinedLaplacian(hObject,handles,1);




function loadniiname2_Callback(hObject, eventdata, handles)
Visualize(hObject,handles,8);



% --- Executes during object creation, after setting all properties.
function loadniiname2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadniiname2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function loadniiname1_Callback(hObject, eventdata, handles)
Visualize(hObject,handles,7);



% --- Executes during object creation, after setting all properties.
function loadniiname1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function QSMProcessing_CreateFcn(hObject, eventdata, handles)

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
if value >SS(3)
    value =SS(3);
end
set(handles.FrameSelect,'value',value);
guidata(hObject, handles);
SliceSelection(hObject,handles);

function savemat1_Callback(hObject, eventdata, handles)
Imag1=handles.Mi;
seed=get(handles.loadniiname1,'string');
uisave('Imag1',seed);

function savemat2_Callback(hObject, eventdata, handles)
Imag1=handles.M2;
seed=get(handles.loadniiname2,'string');
uisave('Imag1',seed);


function loadmat1_Callback(hObject, eventdata, handles)
LoadImages(hObject,handles,4);

function loadmat2_Callback(hObject, eventdata, handles)
LoadImages(hObject,handles,5);



function PermuteImage_Callback(hObject, eventdata, handles)
PermuteImages(hObject,handles,1)

function flipdim_Callback(hObject, eventdata, handles)
PermuteImages(hObject,handles,2)



% --- Executes on button press in PhaseScale2Pi.
function PhaseScale2Pi_Callback(hObject, eventdata, handles)
PermuteImages(hObject,handles,3)




function ZeroPaddingAroundBrain_Callback(hObject, eventdata, handles)
% hObject    handle to ZeroPaddingAroundBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ZeroPaddingAroundBrain as text
%        str2double(get(hObject,'String')) returns contents of ZeroPaddingAroundBrain as a double


% --- Executes during object creation, after setting all properties.
function ZeroPaddingAroundBrain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZeroPaddingAroundBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PhasePaddingsize_Callback(hObject, eventdata, handles)
% hObject    handle to PhasePaddingsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PhasePaddingsize as text
%        str2double(get(hObject,'String')) returns contents of PhasePaddingsize as a double


% --- Executes during object creation, after setting all properties.
function PhasePaddingsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhasePaddingsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getB0Direction.
function getB0Direction_Callback(hObject, eventdata, handles)
 msgbox('Please use axial scan, i.e. B0_dir =[0 0 1]. For oblique scan, please contact us to get the correct method to calculate B0_dir, otherwise the result will be wrong :(')





% --- Executes when selected object is changed in Left.
function Left_SelectionChangeFcn(hObject, eventdata, handles)
UniqueID=get(handles.UniqueID,'string');
try
    switch get (hObject,'Tag')
        case 'L1'
            %handles.Mi=evalin('base', [UniqueID '_magni']);
            handles.Mi=evalin('base', [UniqueID '_magni']);
            set (handles.loadniiname1,'string',[UniqueID '_magni']);
        case 'L2'
            %handles.Mi=evalin('base', [UniqueID '_phase']);
            handles.Mi=evalin('base', get(handles.RightImageLoad,'string'));
            set (handles.loadniiname1,'string',[UniqueID '_phase']);
        case 'L3'
           % handles.Mi=evalin('base', [UniqueID '_Mask']);
            handles.Mi=single(evalin('base', get(handles.NewMaskName,'string')));
            set (handles.loadniiname1,'string',[UniqueID '_Mask']);
        case 'L4'
            handles.Mi=evalin('base', [UniqueID '_Laplacian']);
            set (handles.loadniiname1,'string',[UniqueID '_Laplacian']);
        case 'L5'
            handles.Mi=evalin('base', [UniqueID '_TissuePhase']);
            set (handles.loadniiname1,'string',[UniqueID '_TissuePhase']);
        case 'L6'
            handles.Mi=evalin('base', [UniqueID '_TissuePhase_V']);
            set (handles.loadniiname1,'string',[UniqueID '_TissuePhase_V']);
        case 'L7'
            handles.Mi=evalin('base', [UniqueID '_FrequencyShift']);
             set (handles.loadniiname1,'string',[UniqueID '_FrequencyShift']);
       case 'L8'
            handles.Mi=evalin('base', [UniqueID '_QSM_iLSQR']);
            set (handles.loadniiname1,'string',[UniqueID '_QSM_iLSQR']);
       case 'L9'
            handles.Mi=evalin('base', [UniqueID '_QSM_Fast']);
            set (handles.loadniiname1,'string',[UniqueID '_QSM_Fast']);
    end

    guidata(hObject, handles);
    SliceSelection(hObject,handles);
    AdjustIntensity(hObject,handles,3)
end


% --- Executes when selected object is changed in Right.
function Right_SelectionChangeFcn(hObject, eventdata, handles)
UniqueID=get(handles.UniqueID,'string');
try
    switch get (hObject,'Tag')
        case 'R1'
            handles.M2=evalin('base', [UniqueID '_magni']);
            set (handles.loadniiname2,'string',[UniqueID '_magni']);
        case 'R2'
            handles.M2=evalin('base', [UniqueID '_phase']);
            set (handles.loadniiname2,'string',[UniqueID '_phase']);
       case 'R3'
            handles.M2=evalin('base', [UniqueID '_Mask']);
            set (handles.loadniiname2,'string',[UniqueID '_Mask']);
        case 'R4'
            handles.M2=evalin('base', [UniqueID '_Laplacian']);
            set (handles.loadniiname2,'string',[UniqueID '_Laplacian']);
        case 'R5'
            handles.M2=evalin('base', [UniqueID '_TissuePhase']);
            set (handles.loadniiname2,'string',[UniqueID '_TissuePhase']);
      case 'R6'
            handles.M2=evalin('base', [UniqueID '_TissuePhase_V']);
            set (handles.loadniiname2,'string',[UniqueID '_TissuePhase_V']);
       case 'R7'
            handles.M2=evalin('base', [UniqueID '_FrequencyShift']);
            set (handles.loadniiname2,'string',[UniqueID '_FrequencyShift']);
       case 'R8'
            handles.M2=evalin('base', [UniqueID '_QSM_iLSQR']);
            set (handles.loadniiname2,'string',[UniqueID '_QSM_iLSQR']);
       case 'R9'
            handles.M2=evalin('base', [UniqueID '_QSM_Fast']);
            set (handles.loadniiname2,'string',[UniqueID '_QSM_Fast']);
  end
    guidata(hObject, handles);
    SliceSelection(hObject,handles);
    AdjustIntensity(hObject,handles,3)
end



% --- Executes during object creation, after setting all properties.
function ImageAxes_CreateFcn(hObject, eventdata, handles)



% --- Executes on button press in loadSpatialRes.
function loadSpatialRes_Callback(hObject, eventdata, handles)
LoadImages(hObject,handles,8);




function Voxelsize5_Callback(hObject, eventdata, handles)
% hObject    handle to Voxelsize5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Voxelsize5 as text
%        str2double(get(hObject,'String')) returns contents of Voxelsize5 as a double


% --- Executes during object creation, after setting all properties.
function Voxelsize5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Voxelsize5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
