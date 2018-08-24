function LoadUserData(hObject,handles)
% Wei Li, Duke University, June 2010

[file,file_path] = uigetfile('*.mat');
my_var=load([file_path file],'Susceptibility','R2Star','magni');
Images=my_var.Susceptibility;

fileName=strcat(file_path,file);
[~,f,e] = fileparts(fileName);

         
handles.Mi=permute(Images(:,end:-1:1,:,:),[2,1,3]);
handles.susceptibility=permute(Images(:,end:-1:1,:,:),[2,1,3]);
handles.R2Star=permute(my_var.R2Star(:,end:-1:1,:,:),[2,1,3]);
handles.magni=permute(my_var.magni(:,end:-1:1,:,:),[2,1,3,4]);

I1.img=[];
handles.header=I1;
        
ss=size(handles.Mi);     

set(handles.FrameSelect, 'Max', ss(3),'SliderStep',[1/ss(3) 2/ss(3)]); 
set(handles.FrameSelect, 'visible','on'); 

handles.M2=ones(size(handles.Mi,1),size(handles.Mi,2),size(handles.Mi,3));
handles.img1(:,:,1)=real(handles.M2(:,:,round(ss(3)/2)));
handles.img1(:,:,3)=0;
handles.himage1=imshow(handles.img1,'parent',handles.ImageAxes);
handles.himage2=imshow(real(handles.Mi(:,:,round(ss(3)/2))),[],'parent',handles.ImageAxes2);
AlphaMap=handles.M2(:,:,round(ss(3)/2));
AlphaMap=(AlphaMap>0);
set(handles.himage2,'AlphaData',AlphaMap);
set(handles.himage1,'CDataMapping','scaled')
if get(handles.Grey,'value')
    colormap gray
else 
    colormap jet
end

sss=size(handles.Mi,4);
if sss>1
    set(handles.Coilslider,'Max', sss,'SliderStep',[1/sss 2/sss]); 
    if get(handles.Coilslider, 'value')>sss;
        set(handles.Coilslider, 'value', sss/2); 
        
    end
    set(handles.Coilslider,'value', 1); 
end

Xlim2=get (handles.ImageAxes2,'Xlim');
Ylim2=get (handles.ImageAxes2,'ylim');

XYLimMax=max([Xlim2(2)-Xlim2(1) Ylim2(2)-Ylim2(1)]);

set (handles.ImageAxes2,'Xlim',[0 XYLimMax]+0.5,'Ylim',[0 XYLimMax]+0.5);
set (handles.ImageAxes,'Xlim',[0 XYLimMax]+0.5,'Ylim',[0 XYLimMax]+0.5);
set(handles.FrameSelect, 'value', round(ss(3)/2)); 
set(handles.status,'string',f)
handles.ImageAxisNum=3;
handles.NumberSelect=1;
handles.baseAlpha=0.4;

TempROIvalues=ones(size(handles.Mi));
try
for ii=1:9
    MyROIs=(Atlas==ii);
    TempROIvalues(MyROIs)=ii*0.01+handles.baseAlpha;
end
end
handles.M2=TempROIvalues;
handles.savedata=0;
guidata(hObject, handles);

set(handles.FrameSelect,'value',27);
AdjustIntensity(hObject,handles,20)

SliceSelectionROI(hObject,handles);

