function varargout = STISuite(varargin)
% STISUITE MATLAB code for STISuite.fig
%      STISUITE, by itself, creates a new STISUITE or raises the existing
%      singleton*.
%
%      H = STISUITE returns the handle to a new STISUITE or the handle to
%      the existing singleton*.
%
%      STISUITE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STISUITE.M with the given input arguments.
%
%      STISUITE('Property','Value',...) creates a new STISUITE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before STISuite_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to STISuite_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help STISuite

% Last Modified by GUIDE v2.5 22-Jan-2014 10:51:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @STISuite_OpeningFcn, ...
                   'gui_OutputFcn',  @STISuite_OutputFcn, ...
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


% --- Executes just before STISuite is made visible.
function STISuite_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to STISuite (see VARARGIN)

% Choose default command line output for STISuite
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes STISuite wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = STISuite_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in RGBVisualize.
function RGBVisualize_Callback(hObject, eventdata, handles)
% hObject    handle to RGBVisualize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Color_GUI



% --- Executes on button press in GrayScaleVisualize.
function GrayScaleVisualize_Callback(hObject, eventdata, handles)
% hObject    handle to GrayScaleVisualize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Gray_GUI


% --- Executes on button press in Templates.
function Templates_Callback(hObject, eventdata, handles)
% hObject    handle to Templates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open STI_Templates.m


% --- Executes on button press in Introduction.
function Introduction_Callback(hObject, eventdata, handles)
% hObject    handle to Introduction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open STI_SUITE_Readme.m


% --- Executes on button press in GrayScaleROIEditing.
function GrayScaleROIEditing_Callback(hObject, eventdata, handles)
% hObject    handle to GrayScaleROIEditing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Gray_ROI


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
QSM_GUI


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
open STI_Parfor.m

% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
