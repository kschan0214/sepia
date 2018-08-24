function LoadBinFile(hObject,handles,control)
% Wei Li, Duke University, June 2010

switch control
    case 1
        filters = {'*.bin','BIN-files (*.bin)'};
        seed = 'InputPhaseImage.bin';
        [fn,pn,filterindex] = uigetfile(filters, sprintf('Save Workspace Variables'), seed);
        dimention=[256 256 256];
        if fn==0; return;end
        fid = fopen(fn);
        imag = fread(fid, dimention(1)*dimention(2)*dimention(3), 'float');
        fclose(fid);
        handles.Mi=reshape(imag,dimention);
   case 2
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load NIFTI files'));
        if fn==0; return;end
        try
        I1=readFileNifti(strcat(pn,fn));
        I1=I1.data;
        catch
        I1 = load_untouch_nii(strcat(pn,fn));
        I1=I1.img;
        end
        
        handles.Mi=I1;
   case 22
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load NIFTI files'));
        if fn==0; return;end
        I1=readFileNifti(strcat(pn,fn));
        
        str1=sprintf('[ %g  %g  %g ]',I1.pixdim(1:3));
        I1=I1.data;

        set(handles.LeftResolution,'string',str1); 
        set(handles.VoxelSizeForQSM,'string',str1); 
        set(handles.Voxelsize5,'string',str1); 
        
        handles.Mi=I1;
        UniqueID=get(handles.UniqueID,'string');
        set(handles.VariableSelect,'string',[UniqueID '_magni']); 
        assignin('base',[UniqueID '_magni'],handles.Mi);

   case 3
        filters = {'*.nii','NII-files (*.nii)';'*.*','all files(*.*)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('Save Workspace Variables'));
        if fn==0; return;end
        I1=load_nii(strcat(pn,fn));
        I1=I1.img;
        handles.Mi=I1;    
   case 4
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load NIFTI files'));
        if fn==0; return;end
        I1=readFileNifti(strcat(pn,fn));
        I1=I1.data;
        handles.M2=I1;
   case 44
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load NIFTI files'));
        if fn==0; return;end
        I1=readFileNifti(strcat(pn,fn));
        I1=I1.data;
        handles.M2=I1;
        UniqueID=get(handles.UniqueID,'string');
        set(handles.RightImageLoad,'string',[UniqueID '_phase']); 
        assignin('base',[UniqueID '_phase'],handles.M2);
 case 45
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load NIFTI files'));
        if fn==0; return;end
        I1=readFileNifti(strcat(pn,fn));
        I1=double(I1.data);
        
        handles.M2=I1;
        UniqueID=get(handles.UniqueID,'string');
        set(handles.BrainMaskName,'string',[UniqueID '_Mask']); 
        set(handles.NewMaskName,'string',[UniqueID '_Mask']);
        assignin('base',[UniqueID '_Mask'],handles.M2);

   case 5
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('Save Workspace Variables'));
        if fn==0; return;end
        fileName=strcat(pn,fn);
        [p,f,e] = fileparts(fileName);
        %--------------
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
        %-------------
        handles.Mi=I1.img;
        I1.img=[];
        handles.header=I1;
    case 6
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('Save Workspace Variables'));
        if fn==0; return;end
        fileName=strcat(pn,fn);
        [p,f,e] = fileparts(fileName);
        %--------------
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
        %-------------
        handles.M2=I1.img;
        handles.Phase=I1.img;
        I1.img=[];
        handles.header2=I1;
     case 7
        filters = {'*.spr;*.sdt'};
        [fn,pn,~] = uigetfile(filters, sprintf('load sdt files'));
        if fn==0; return;end
        %I1=readFileNifti(strcat(pn,fn));
        filename=strcat(pn,fn);
        
        [I1,~,~,~,~,~]=readsdt(filename(1:end-4));
        I1=permute(I1,[2 1 3 4]);
        I1=I1(end:-1:1,:,:,:);
        handles.Mi=I1;
       %set(handles.VariableSelect,'string',fn);
     case 8
        filters = {'*.spr;*.sdt'};
        [fn,pn,~] = uigetfile(filters, sprintf('load sdt files'));
        if fn==0; return;end
        %I1=readFileNifti(strcat(pn,fn));
        filename=strcat(pn,fn);
        
        [I1,~,~,~,~,~]=readsdt(filename(1:end-4));
        I1=permute(I1,[2 1 3 4]);
        I1=I1(end:-1:1,:,:,:);
        handles.M2=I1;
       %set(handles.VariableSelect,'string',fn);
   case 23
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load NIFTI files'));
        if fn==0; return;end
        I1=readFileNifti(strcat(pn,fn));
        I1=I1.data;
        handles.Mi=I1;
        set(handles.LoadRealData,'string','***'); 

   case 26
        filters = {'*.nii.gz;*.nii','NIFTI-files(*.nii.gz,*.nii)'};
        [fn,pn,filterindex] = uigetfile(filters, sprintf('load NIFTI files'));
        if fn==0; return;end
        I1=readFileNifti(strcat(pn,fn));
        I1=I1.data;
        Imagdata=I1;
        Complexdata=handles.Mi+1i*Imagdata;
        Magdata=abs(Complexdata);
        Phasedata=angle(Complexdata);
        handles.M2=Phasedata;
        handles.Mi=Magdata;
        UniqueID=get(handles.UniqueID,'string');
        assignin('base',[UniqueID '_magni'],Magdata);
        set(handles.VariableSelect,'string',[UniqueID '_magni']); 
        assignin('base',[UniqueID '_phase'],Phasedata);
        set(handles.RightImageLoad,'string',[UniqueID '_phase']); 
        ss=size(handles.Mi);
        Image=real(handles.Mi(:,:,round(ss(3)/2)));
        Image=permute(Image(:,end:-1:1),[2 1]);
        handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
        set(handles.LoadImagData,'string','***'); 
end

if control==4||control==6||control==8||control==44||control==45||control==26
ss=size(handles.Mi);
Image=real(handles.M2(:,:,round(ss(3)/2)));
Image=permute(Image(:,end:-1:1),[2 1]);
handles.himage2=imshow(Image,[],'parent',handles.ImageAxes2);
else
ss=size(handles.Mi);
set(handles.FrameSelect, 'Max', ss(3),'SliderStep',[1/ss(3) 2/ss(3)]); 
Image=real(handles.Mi(:,:,round(ss(3)/2)));
Image=permute(Image(:,end:-1:1),[2 1]);
handles.himage1=imshow(Image,[],'parent',handles.ImageAxes);
set(handles.FrameSelect, 'value', round(ss(3)/2)); 
end
set(handles.status, 'String', [pn fn]); 
sss=size(handles.Mi,4);
if sss>1
    set(handles.Coilslider,'Min',1,'Max', sss,'SliderStep',[1/(sss-1) 2/(sss-1)]); 
    if get(handles.Coilslider, 'value')>sss;
        set(handles.Coilslider, 'value', sss/2); 
    end
    set(handles.Coilslider,'value', 1); 
end

handles.ImageAxisNum=3;
guidata(hObject, handles);
