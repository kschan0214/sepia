function LoadBinFileROI(hObject,handles,control)
% Wei Li, Duke University, June 2010
%control=0;
handles.baseAlpha=0.4;

switch control
    case 0
        filters = {'*.nii.gz;*.nii;*ROI.mat','NIFTI-files(*.nii.gz,*.nii,*ROI.mat)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('Save Workspace Variables'));
        if fn==0; return;end
        
        fileName1=strcat(pn,[fn(1:end-7) 'ROIs_warp']);     
        ROIfile=load (fileName1);
        ROIs_warp=ROIfile.ROIs_warp;
        
        fileName=strcat(pn,fn);
        [~,f,e] = fileparts(fileName);
        if(strcmpi(e,'.gz'))
           tmpDir = tempname;
           mkdir(tmpDir);
           gunzip(strcat(pn,fn), tmpDir);
           tmpFileName = fullfile(tmpDir, f);
           tmpFile = true;
           I1 = load_untouch_nii(tmpFileName);
        else
           tmpFile = false;
           I1 = load_untouch_nii(strcat(pn,fn));
        end
        if(tmpFile)
            delete(tmpFileName);
        end
        Images=I1.img;
    case 1
        [file,file_path,FilterIndex] = uigetfile('*.mat','MultiSelect','on');

        try
            for ii=1:length(file)
                my_var=load([file_path file{ii}]);
                my_name=fieldnames(my_var);
                Images(:,:,:,ii)=getfield(my_var,my_name{1});
            end
        catch
            my_var=load([file_path file]);
            my_name=fieldnames(my_var);
            Images=getfield(my_var,my_name{1});
        end
        ROIs_warp=Images*0;
        f='Images';
        
        
     case 2
         filters = {'*.spr;*.sdt'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load sdt files'));
        if fn==0; return;end
        filename=strcat(pn,fn);
        [I1,NS,NP,NI,NR,precision]=readsdt(filename(1:end-4));
        Images=I1;
        ROIs_warp=Images*0;
        f=fn;
     case 5
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [file,file_path] = uigetfile(filters, sprintf('load NIFTI files'),'MultiSelect','on');
        try
            for ii=1:length(file)
                I1=readFileNifti(strcat(file_path,file{ii}));
                Images(:,:,:,ii)=I1.data;
            end
        catch
            I1=readFileNifti(strcat(file_path,file));
            Images=I1.data;
        end
        ROIs_warp=Images*0;
        f='Images';
        
end

handles.Mi=Images;
ROIs_warp1=ROIs_warp;

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
    set(handles.Coilslider,'Max', sss,'SliderStep',[1/(sss-1) 2/(sss-1)]); 
    if get(handles.Coilslider, 'value')>sss;
        set(handles.Coilslider, 'value', round(sss/2)); 
    end
    set(handles.Coilslider,'value', 1); 
end

Xlim2=get (handles.ImageAxes2,'Xlim');
Ylim2=get (handles.ImageAxes2,'ylim');

XYLimMax=max([Xlim2(2)-Xlim2(1) Ylim2(2)-Ylim2(1)]);

set (handles.ImageAxes2,'Xlim',[0 XYLimMax]+0.5,'Ylim',[0 XYLimMax]+0.5);
set (handles.ImageAxes,'Xlim',[0 XYLimMax]+0.5,'Ylim',[0 XYLimMax]+0.5);
set(handles.FrameSelect, 'value', round(ss(3)/2)); 
set (handles.ImageAxes2,'ydir','normal');
set (handles.ImageAxes,'ydir','normal');

set(handles.status,'string',f)
handles.ImageAxisNum=3;
handles.NumberSelect=1;

guidata(hObject, handles);

%set(handles.FrameSelect,'value',27);
SliceSelectionROI(hObject,handles);

TempROIvalues=ones(size(handles.M2));
for ii=1:9
    MyROIs=(ROIs_warp1==ii);
    TempROIvalues(MyROIs)=ii*0.01+handles.baseAlpha;
end            
handles.M2=TempROIvalues;
guidata(hObject, handles);
        
