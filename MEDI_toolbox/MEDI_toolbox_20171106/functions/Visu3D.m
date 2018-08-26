function varargout = Visu3D(varargin)
% VISU3D M-file for Visu3D.fig
%      VISU3D, by itself, creates a new VISU3D or raises the existing
%      singleton*.
%
%      H = VISU3D returns the handle to a new VISU3D or the handle to
%      the existing singleton*.
%
%      VISU3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISU3D.M with the given input arguments.
%
%      VISU3D('Property','Value',...) creates a new VISU3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Visu3D_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Visu3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help Visu3D

% Last Modified by GUIDE v2.5 12-Dec-2007 12:20:52

% Begin initialization code - DO NOT EDIT

%
%   example: Visu3D(QSM,'dimension',voxel_size,'defaultposition',[100 100 100])
% 
%   Created by Ludovic de Rochefort
%   Last modified by Pascal Spincemaille

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Visu3D_OpeningFcn, ...
                   'gui_OutputFcn',  @Visu3D_OutputFcn, ...
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


% --- Executes just before Visu3D is made visible.
function Visu3D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Visu3D (see VARARGIN)

if size(varargin,2)>0
handles.iM=varargin{1};
handles.iM=handles.iM(:,:,:,:);
else
    load mri
    handles.iM=double(permute(D,[2 1 4 3]));
end

handles.dimension=[1 1 2.5];
handles.Show_pos=0;
handles.Wait=0;
handles.save_im=0;
handles.filen='';

nam=fieldnames(handles);
if sum(strcmp(nam,'x'))==0
handles.x=floor((size(handles.iM,1)+1)/2);
handles.y=floor((size(handles.iM,2)+1)/2);
handles.z=floor((size(handles.iM,3)+1)/2);
end

if handles.x>size(handles.iM,1)
handles.x=floor((size(handles.iM,1)+1)/2);
handles.y=floor((size(handles.iM,2)+1)/2);
handles.z=floor((size(handles.iM,3)+1)/2);    
end


if size(varargin,2)>1
    for k=1:(size(varargin,2)-1)/2
        if strcmpi(varargin{2*k},'dimension')
            handles.dimension=varargin{2*k+1};
        end
        if strcmpi(varargin{2*k},'showpos')
            handles.Show_pos=varargin{2*k+1};
            set(handles.show_pos,'Value',varargin{2*k+1});
        end
        if strcmpi(varargin{2*k},'absolute')
            set(handles.abs,'Value',varargin{2*k+1});
            if varargin{2*k+1}
            set(handles.plot_type,'Value',3);
            end
        end
        if strcmpi(varargin{2*k},'colormap')
            set(handles.colormap,'Value',varargin{2*k+1});
        end
        if strcmpi(varargin{2*k},'Wait')
            handles.Wait=varargin{2*k+1};
        end
        if strcmpi(varargin{2*k},'DefaultPosition')
            pos=varargin{2*k+1};
            handles.x=min(max(pos(1),1),size(handles.iM,1));
            handles.y=min(max(pos(2),1),size(handles.iM,2));
            handles.z=min(max(pos(3),1),size(handles.iM,3));
        end
    end
end


handles.Position=[handles.x handles.y handles.z];

if handles.Wait
    set(handles.figure1,'CloseRequestFcn','uiresume(gcbo); set(gcbo,''CloseRequestFcn'',''closereq'')')
end
set(handles.figure1, 'DoubleBuffer', 'on');

set(handles.coronal,'visible','off')
set(handles.sagittal,'visible','off')
set(handles.axial,'visible','off')
set(handles.colorbar,'visible','off')

set(handles.echo,'max',size(handles.iM,4)+0.5);
set(handles.echo,'value',1);
set(handles.echo_edit,'string','1');

if size(handles.iM,4)>1
    set(handles.echo,'sliderstep',[1/(size(handles.iM,4)-1) 1/(size(handles.iM,4)-1)]);
end
handles.echo_num=1;

handles.zoom=1;
handles.min=0;
handles.max=1;

ab=abs(handles.iM(:));
handles.ref_max=max(ab(~isnan(ab)));
if size(handles.ref_max,1)==0
    handles.ref_max=1;
end

handles.ab2=[max(real(handles.iM(:))) max(imag(handles.iM(:))) ...
    max(abs(handles.iM(:))) max(angle(handles.iM(:)))];

% This makes the scaling slider produce images in a new figure instead of
% the GUI's figure. Need to figure out why
% hListener=handle.listener(handles.scaling_slider,'ActionEvent', @scaling_slider_Callback);
% setappdata(handles.scaling_slider,'myListener',hListener);
% Choose default command line output for Visu3D
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

plot_image(hObject, eventdata, handles);

% UIWAIT makes Visu3D wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = Visu3D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout(1) = {[handles.x handles.y handles.z]};
varargout{1} = handles.output;


% --- Executes on button press in autolevel.
function autolevel_Callback(hObject, eventdata, handles)
% hObject    handle to autolevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autolevel

handles.Autolevel=get(hObject,'Value');



% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in show_pos.
function show_pos_Callback(hObject, eventdata, handles)
% hObject    handle to show_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_pos

handles.Show_pos=get(hObject,'Value');
plot_image(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);





% --- Executes on slider movement.
function echo_Callback(hObject, eventdata, handles)
% hObject    handle to echo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


handles.echo_num=floor(get(hObject,'Value'));
set(handles.echo_edit,'String',num2str(handles.echo_num))
plot_image(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function echo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to echo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function echo_edit_Callback(hObject, eventdata, handles)
% hObject    handle to echo_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of echo_edit as text
%        str2double(get(hObject,'String')) returns contents of echo_edit as a double



handles.echo_num=floor(str2double(get(hObject,'String')));
handles.echo_num=max(handles.echo_num,get(handles.echo,'min'));
handles.echo_num=floor(min(handles.echo_num,get(handles.echo,'max')));


set(handles.echo_edit,'String',num2str(handles.echo_num))
set(handles.echo,'Value',(handles.echo_num))
plot_image(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function echo_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to echo_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in abs.
function abs_Callback(hObject, eventdata, handles)
% hObject    handle to abs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of abs

plot_image(hObject, eventdata, handles)

% --- Executes on selection change in colormap.
function colormap_Callback(hObject, eventdata, handles)
% hObject    handle to colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns colormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from colormap

contents = get(hObject,'String');
colormap(contents{get(hObject,'Value')});
if get(handles.nan,'value')
    set_up_nan_colormap;
end
% set(handles.figure1,'currentaxes',handles.colorbar);
% colorbar(handles.colorbar,'peer',handles.coronal)



% --- Executes during object creation, after setting all properties.
function colormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_up_nan_colormap

map=get(gcf,'colormap');
map(1,:)=[1 1 1];
set(gcf,'colormap',map)


function [imcor imsag imaxi ma iMsize xyzval] = compute_image(hObject, eventdata, handles)
typ=get(handles.plot_type,'value');
nor=get(handles.norm_type,'value');
f={@real, @imag, @abs, @angle};

iMsize=size(handles.iM);
if nor==1
    ab=handles.ab2(typ);
end
if nor==2
    ab=abs(f{typ}(handles.iM(:,:,:,handles.echo_num)));
end
if nor==3
    abx=abs(f{typ}(handles.iM(handles.x,:,:,handles.echo_num)));
    aby=abs(f{typ}(handles.iM(:,handles.y,:,handles.echo_num)));
    abz=abs(f{typ}(handles.iM(:,:,handles.z,handles.echo_num)));
    ab=[abx(:);aby(:);abz(:)];
end
ma=max(ab(~isnan(ab)));
if size(ma,1)==0
   ma=1;
end
imcor=f{typ}(handles.iM(:,handles.y,end:-1:1,handles.echo_num));
imcor=transpose(squeeze(imcor));

imsag=f{typ}(handles.iM(handles.x,:,end:-1:1,handles.echo_num));
imsag=transpose(squeeze(imsag));

imaxi=f{typ}(handles.iM(:,:,handles.z,handles.echo_num));
imaxi=transpose(squeeze(imaxi));
imaxi(end:-1:1,:)=imaxi;
xyzval=handles.iM(handles.x,handles.y,handles.z,handles.echo_num);

function plot_image(hObject, eventdata, handles)
% handles=guidata(hObject);

% handles.iM=1e6*handles.iM;
scale=get(handles.scaling_slider,'value');
win_level=get(handles.win_level_slider,'value');

% [scale win_level]

if scale==0
   scale=eps; 
end

save_im=handles.save_im;
filen=handles.filen;
[path,filen,ext]=fileparts(filen);
filen=[filen '_' num2str(handles.x) '_' num2str(handles.y) '_' num2str(handles.z)];
   
contents = get(handles.colormap,'String');

% if get(handles.abs,'Value')
%     dr=([0 (handles.max)*handles.ref_max]);
%     if dr(2)==0;
%         dr(2)=1;
%     end
% else
%     dr=[-(handles.max)*(handles.ref_max) (handles.max)*handles.ref_max];
%     if dr(2)==0;
%         dr=[-1 1];
%     end
% end


[imcor imsag imaxi ma iMsize xyzval] = compute_image(hObject, eventdata, handles);

cont=get(handles.contour,'value');
lev=eval(get(handles.contour_level,'string'));

txt={['Location: ' num2str([handles.x handles.y handles.z])]};
txt(1)={[txt{1} ', Value: ' num2str(xyzval)]};
set(handles.stat,'value',1)
set(handles.stat,'String',txt)

% hold off



axx=[-ma ma];
if get(handles.abs,'Value')
    axx=[0 ma];
end

axx=axx*scale;
axx=axx+ma*win_level;

% iM(iM==0)=NaN;

axes(handles.coronal);

if sum(isnan(imcor(:)))~=prod(size(imcor))
    % imshow(im,'Displayrange',dr)
    imshow(imcor)
end
caxis(axx)
daspect([handles.dimension(3) handles.dimension(1) 1]/max(handles.dimension))
% hold on
if handles.Show_pos
    plot([0.5 iMsize(1)+0.5],iMsize(3)+0.5-[handles.z handles.z],'-r')
    plot([handles.x handles.x],[0.5 iMsize(3)+0.5],'-r')
end

if cont
    [X,Y]=meshgrid(1:size(imcor,2),1:size(imcor,1));
    v=floor((min(imcor(:)):lev:lev+max(imcor(:)))/lev)*lev;
    if (size(imcor,1)>2).*(size(imcor,2)>2)
    contour(X,Y,imcor,v,'-k');
    end
end


if save_im
    tag='coronal'
    h=figure('name',tag)
    % imshow(imcor,'Displayrange',dr)
    imshow(imcor)
    caxis(axx)
    daspect(gca,[handles.dimension(3) handles.dimension(1) 1]/max(handles.dimension))
    hold on
    if handles.Show_pos
        plot([0.5 iMsize(1)+0.5],iMsize(3)+0.5-[handles.z handles.z],'-r')
        plot([handles.x handles.x],[0.5 iMsize(3)+0.5],'-r')
    end
    
    if cont
        [X,Y]=meshgrid(1:size(imcor,2),1:size(imcor,1));
        v=floor((min(imcor(:)):lev:lev+max(imcor(:)))/lev)*lev;
        if (size(imcor,1)>2).*(size(imcor,2)>2)
            contour(X,Y,imcor,v,'-k');
        end
    end
    
    set(gcf,'color','white')
    a=0.05;
    % [path,filen,ext]=fileparts(filen);
    % tag='visu3d'
    set(gcf,'units','normalized')
    set(gcf,'position',[a a 1-2*a 1-2*a])
    set(gcf,'units','pixel')
    % if strcmp(path,'')==1
    %     path='.';
    % end
    %     pos=get(gcf,'position')
    %     set(gcf,'position',[pos([1 2]) size(imcor,2) size(imcor,1)])
    set(gca,'position',[0 0 1 1])
    colormap(contents{get(handles.colormap,'Value')});
    if get(handles.nan,'value')
        set_up_nan_colormap;
    end
    drawnow
    saveas(h,[path '/' filen '_' tag ext])
    close(h)
end



axes(handles.sagittal);

if sum(isnan(imsag(:)))~=prod(size(imsag))
% imshow(imsag,'Displayrange',dr)
imshow(imsag)

end

caxis(axx)
daspect([handles.dimension(3) handles.dimension(2) 1]/max(handles.dimension))
% hold on
if handles.Show_pos
    plot([0.5 iMsize(2)+0.5],iMsize(3)+0.5-[handles.z handles.z],'-r')
    plot([handles.y handles.y],[0.5 iMsize(3)+0.5],'-r')
end
if cont
    [X,Y]=meshgrid(1:size(imsag,2),1:size(imsag,1));
    v=floor((min(imsag(:)):lev:lev+max(imsag(:)))/lev)*lev;
    if (size(imsag,1)>2).*(size(imsag,2)>2)
        contour(X,Y,imsag,v,'-k');
    end
end

if save_im
    tag='sagittal'
    h=figure('name',tag)
    % imshow(imsag,'Displayrange',dr)
    imshow(imsag)
    
    caxis(axx)
    daspect([handles.dimension(3) handles.dimension(2) 1]/max(handles.dimension))
    hold on
    if handles.Show_pos
        plot([0.5 iMsize(2)+0.5],iMsize(3)+0.5-[handles.z handles.z],'-r')
        plot([handles.y handles.y],[0.5 iMsize(3)+0.5],'-r')
    end
    
    if cont
        [X,Y]=meshgrid(1:size(imsag,2),1:size(imsag,1));
        v=floor((min(imsag(:)):lev:lev+max(imsag(:)))/lev)*lev;
        if (size(imsag,1)>2).*(size(imsag,2)>2)
            contour(X,Y,imsag,v,'-k');
        end
    end

    set(h,'color','white')
    a=0.05;
    % [path,filen,ext]=fileparts(filen);
    % tag='visu3d'
    set(h,'units','normalized')
    set(h,'position',[a a 1-2*a 1-2*a])
    set(h,'units','pixel')
    % if strcmp(path,'')==1
    %     path='.';
    % end
    %     pos=get(gcf,'position')
    %     set(gcf,'position',[pos([1 2]) size(imsag,2) size(imsag,1)])
    set(gca,'position',[0 0 1 1])
    colormap(contents{get(handles.colormap,'Value')});
    if get(handles.nan,'value')
        set_up_nan_colormap;
    end
    drawnow
    saveas(h,[path '/' filen '_' tag ext])
    close(h)
    
end

axes(handles.axial);
if sum(isnan(imaxi(:)))~=prod(size(imaxi))
    % imshow(imaxi,'Displayrange',dr)
    imshow(imaxi)
    
end

caxis(axx)
daspect([handles.dimension(2) handles.dimension(1) 1]/max(handles.dimension))
% hold on
if handles.Show_pos
    plot([0.5 iMsize(1)+0.5],iMsize(2)+1-[handles.y handles.y],'-r')
    plot([handles.x handles.x],[0.5 iMsize(2)+0.5],'-r')
end

if cont
    [X,Y]=meshgrid(1:size(imaxi,2),1:size(imaxi,1));
    v=floor((min(imaxi(:)):lev:lev+max(imaxi(:)))/lev)*lev;
    if (size(imaxi,1)>2).*(size(imaxi,2)>2)
        contour(X,Y,imaxi,v,'-k');
    end
end

if save_im
    tag='axial'
    h=figure('name',tag)
    % imshow(imaxi,'Displayrange',dr)
    hh=imshow(imaxi);
    
    caxis(axx)
    daspect([handles.dimension(2) handles.dimension(1) 1]/max(handles.dimension))
    hold on
    if handles.Show_pos
        plot([0.5 iMsize(1)+0.5],iMsize(2)+1-[handles.y handles.y],'-r')
        plot([handles.x handles.x],[0.5 iMsize(2)+0.5],'-r')
    end
    
    if cont
        [X,Y]=meshgrid(1:size(imaxi,2),1:size(imaxi,1));
        v=floor((min(imaxi(:)):lev:lev+max(imaxi(:)))/lev)*lev;
        contour(X,Y,imaxi,v,'-k');
    end
    
    set(h,'color','white')
    a=0.05;
    % [path,filen,ext]=fileparts(filen);
    % tag='visu3d'
    set(h,'units','normalized')
    set(h,'position',[a a 1-2*a 1-2*a])
    set(h,'units','pixel')
    % if strcmp(path,'')==1
    %     path='.';
    % end
    %     pos=get(gcf,'position')
    %     set(gcf,'position',[pos([1 2]) size(imaxi,2) size(imaxi,1)])
    set(gca,'position',[0 0 1 1])
    colormap(contents{get(handles.colormap,'Value')});
    if get(handles.nan,'value')
        set_up_nan_colormap;
    end
    drawnow
    saveas(h,[path '/' filen '_' tag ext])
    close(h)
    
end

% plot colorbar
set(handles.figure1,'currentaxes',handles.colorbar);

% caxis(dr);
% hold on
% imshow(0,'Displayrange',dr)
caxis(axx)
colorbar;


% colorbar(handles.colorbar,'peer',handles.coronal)
% colorbar(handles.colorbar)
contents = get(handles.colormap,'String');
colormap(contents{get(handles.colormap,'Value')});
if get(handles.nan,'value')
    set_up_nan_colormap;
end
% set(handles.figure1,'units','pixel')
% pixval on

if save_im
    tag='cmap'
    h=figure('name',tag)
    caxis(axx)
    colorbar;


    % colorbar(handles.colorbar,'peer',handles.coronal)
    % colorbar(handles.colorbar)
    contents = get(handles.colormap,'String');
    colormap(contents{get(handles.colormap,'Value')});
    if get(handles.nan,'value')
        set_up_nan_colormap;
    end
    
    set(gcf,'color','white')
    a=0.05;
    % [path,filen,ext]=fileparts(filen);
    % tag='visu3d'
    set(h,'units','normalized')
    set(h,'position',[a a 1-2*a 1-2*a])
    set(h,'units','pixel')
    % if strcmp(path,'')==1
    %     path='.';
    % end
    drawnow
    saveas(h,[path '/' filen '_' tag ext])
    close(h)
end



% get(handles.coronal,'ButtonDownFcn')
% set(handles.coronal,'ButtonDownFcn',{'Button_coronal(handles)'});
% axes(handles.coronal)

% function Button_coronal(handles)
% 
% disp('OK')


% add plot of every image in one single figure



% --- Executes on button press in move.
function move_Callback(hObject, eventdata, handles)
% hObject    handle to move (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pt=ginput(1);
p=get(handles.figure1,'currentaxes');
c=(p==handles.coronal)+(p==handles.sagittal)+(p==handles.axial);

while c
if p==handles.coronal
    handles.x=floor(pt(1)+0.5);
    handles.z=size(handles.iM,3)+1-floor(pt(2)+0.5);
    plot_image(hObject, eventdata, handles)
end
if p==handles.sagittal
    handles.y=floor(pt(1)+0.5);
    handles.z=size(handles.iM,3)+1-floor(pt(2)+0.5);
    plot_image(hObject, eventdata, handles)
end
if p==handles.axial
    handles.x=floor(pt(1)+0.5);
    handles.y=size(handles.iM,2)+1-floor(pt(2)+0.5);
    plot_image(hObject, eventdata, handles)
end

pt=ginput(1);
p=get(handles.figure1,'currentaxes');
c=(p==handles.coronal)+(p==handles.sagittal)+(p==handles.axial);

end

handles.Position=[handles.x handles.y handles.z];

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'CloseRequestFcn','closereq')

if handles.Wait
    uiresume(handles.figure1);
else
    close(handles.figure1);
end



% --- Executes on button press in contour.
function contour_Callback(hObject, eventdata, handles)
% hObject    handle to contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of contour

plot_image(hObject, eventdata, handles)

function contour_level_Callback(hObject, eventdata, handles)
% hObject    handle to contour_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contour_level as text
%        str2double(get(hObject,'String')) returns contents of contour_level as a double

if get(handles.contour,'value')
plot_image(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function contour_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contour_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on slider movement.
function scaling_slider_Callback(hObject, eventdata, handles)
% hObject    handle to scaling_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

 handles=guidata(hObject);
val=get(handles.scaling_slider,'value');
set(handles.scaling_factor,'string',num2str(val));

% Update handles structure
guidata(hObject, handles);

plot_image(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function scaling_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaling_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function scaling_factor_Callback(hObject, eventdata, handles)
% hObject    handle to scaling_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scaling_factor as text
%        str2double(get(hObject,'String')) returns contents of scaling_factor as a double

val=eval(get(handles.scaling_factor,'string'));
val=min([val,1]);
val=max([val,0]);
set(handles.scaling_factor,'string',num2str(val));
set(handles.scaling_slider,'value',val);

% Update handles structure
guidata(hObject, handles);

plot_image(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function scaling_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaling_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[f,p]=uiputfile('*.fig');

saveas(handles.figure1,[p f]);



% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[f,p]=uigetfile('*.fig');
uiopen([p f],1);
close(handles.figure1);




% --- Executes on button press in save_image.
function save_image_Callback(hObject, eventdata, handles)
% hObject    handle to save_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fo=imformats;
fo=[fo.ext];
for k=1:size(fo,2)
   fo{k}=['*.' fo{k}]; 
end

[f,p]=uiputfile(fo');
if isnumeric(f)==0
handles.save_im=1;
handles.filen=[p,f];
plot_image(hObject, eventdata, handles);
handles.save_im=0;

% Update handles structure
guidata(hObject, handles);
end



% --- Executes on slider movement.
function win_level_slider_Callback(hObject, eventdata, handles)
% hObject    handle to win_level_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=get(handles.win_level_slider,'value');
set(handles.win_level,'string',num2str(val));
set(handles.win_level_slider,'value',val);

% Update handles structure
guidata(hObject, handles);

plot_image(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function win_level_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to win_level_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function win_level_Callback(hObject, eventdata, handles)
% hObject    handle to win_level_slider_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of win_level_slider_level as text
%        str2double(get(hObject,'String')) returns contents of win_level_slider_level as a double
val=eval(get(handles.win_level,'string'));
val=min([val,1]);
val=max([val,-1]);
set(handles.win_level,'string',num2str(val));
set(handles.win_level_slider,'value',val);

% Update handles structure
guidata(hObject, handles);

plot_image(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function win_level_CreateFcn(hObject, eventdata, handles)
% hObject    handle to win_level_slider_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in roi.
function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ca=gca;


typ=get(handles.plot_type,'value');
if typ==1
    iM=real(handles.iM);
end
if typ==2
    iM=imag(handles.iM);
end
if typ==3
    iM=abs(handles.iM);
end
if typ==4
    iM=angle(handles.iM);
end

iM=iM(:,:,:,handles.echo_num);


ok=0;
if ca==handles.coronal
im=iM(:,handles.y,end:-1:1,handles.echo_num);
im=transpose(squeeze(im));
ok=1;
end

if ca==handles.sagittal
im=iM(handles.x,:,end:-1:1,handles.echo_num);

%im=iM(handles.x,end:-1:1,end:-1:1,handles.echo_num); %Tian Liu 2008.02.19

im=transpose(squeeze(im));
ok=1;
end

if ca==handles.axial
im=iM(:,:,handles.z,handles.echo_num);
im=transpose(squeeze(im));
im(end:-1:1,:)=im;
ok=1;
end

if ok==0
    axes(handles.coronal);
im=iM(:,handles.y,end:-1:1,handles.echo_num);
im=transpose(squeeze(im));
    ok=1;
end
    
    
if ok
% [x,y,a,hasimage] = getimage;
[xi,yi] = getline(gcf,'closed'); 
bw=roipoly(im,xi,yi);
plot(xi,yi)


dat=im(bw);
txt=get(handles.stat,'string');
k=size(txt,1);
pos=mean([xi yi]);
text(pos(1),pos(2),num2str(k))
txt(end+1)={['ROI ' num2str(size(txt,1)) ': Mean: ' num2str(mean(dat)) ' Std: ' num2str(std(dat)) ' Size: ' num2str(sum(bw(:)))]};
set(handles.stat,'String',txt)
end
% imshow(bw)
% 
% ROI(im)




% --- Executes on selection change in norm_type.
function norm_type_Callback(hObject, eventdata, handles)
% hObject    handle to norm_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns norm_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from norm_type


    plot_image(hObject, eventdata, handles)

% --- Executes on selection change in plot_type.
function plot_type_Callback(hObject, eventdata, handles)
% hObject    handle to plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns plot_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_type

    plot_image(hObject, eventdata, handles)
    
    
% --- Executes during object creation, after setting all properties.
function plot_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in reset_win_level.
function reset_win_level_Callback(hObject, eventdata, handles)
% hObject    handle to reset_win_level (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val=0;
val=min([val,1]);
val=max([val,0]);
set(handles.win_level,'string',num2str(val));
set(handles.win_level_slider,'value',val);

% Update handles structure
guidata(hObject, handles);

plot_image(hObject, eventdata, handles);

% --- Executes on button press in reset_scaling.
function reset_scaling_Callback(hObject, eventdata, handles)
% hObject    handle to reset_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val=1;
val=min([val,1]);
val=max([val,0]);
set(handles.scaling_factor,'string',num2str(val));
set(handles.scaling_slider,'value',val);

% Update handles structure
guidata(hObject, handles);

plot_image(hObject, eventdata, handles);



% --- Executes on selection change in roi_list.
function roi_list_Callback(hObject, eventdata, handles)
% hObject    handle to roi_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns roi_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from roi_list


% --- Executes during object creation, after setting all properties.
function roi_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in nan.
function nan_Callback(hObject, eventdata, handles)
% hObject    handle to nan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nan

plot_image(hObject, eventdata, handles)







% --- Executes on button press in invert.
function invert_Callback(hObject, eventdata, handles)
% hObject    handle to invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.iM=-handles.iM;
guidata(hObject, handles);
plot_image(hObject, eventdata, handles)