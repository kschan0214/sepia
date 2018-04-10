%%%% Copyright ?WANGLAB 2011 All Rights Reserved. 
% Created by Shuai on 2011.06.13
% for reading bruker data multi-echo 3D data in one file.
% if the multi-echo data is in different file,please use
% Read_Bruker_3D_MultiFile 
% Modified by Alexey on 12/10/2015

function [iField voxel_size matrix_size TE delta_TE CF Affine3D B0_dir TR NumEcho] = Read_Bruker_raw(fid_dir)

fid = fopen([fid_dir '/acqp'],'r');
f1 = fread(fid,Inf,'*char');
fclose(fid);
f1 =f1';

fid = fopen([fid_dir '/method'],'r');
f2 = fread(fid,Inf,'*char');
fclose(fid);
f2 =f2';


ind1 = strfind(f2,'##$PVM_SpatResol');
if (isempty(ind1))
    fprintf('sorry,Can not find "%s" Value,Please check Method file.\n','PVM_SpatResol');
    return;
end
p = f2(ind1+4:end);
ind1 = strfind(p,')');
i1 = ind1(1)+1;
ind1 = strfind(p,'##$');
i2 = ind1(1)-1;
voxel_size = str2num(p(i1:i2));


ind1 = strfind(f2,'##$PVM_Matrix');
if (isempty(ind1))
    fprintf('sorry,Can not find "%s" Value,Please check Method file.\n','PVM_Matrix');
    return;
end
p = f2(ind1+4:end);
ind1 = strfind(p,')');
i1 = ind1(1)+1;
ind1 = strfind(p,'##$');
ind12 = strfind(p,'$$');
i2 = min(ind1(1),ind12(1))-1;
matrix_size = str2num(p(i1:i2));

%%
%multi echo
%##$EffectiveTE=( 12 )
%7.9 15.8 23.7 31.6 39.5 47.4 55.3 63.2 71.1 79 86.9 94.8

ind0 = strfind(f2,'##$PVM_SPackArrPhase1Offset');
p = f2(ind0+4:end);
ind0 = strfind(p,'$$');
ind0 = ind0(1)-1;
ind02 = strfind(p,')');
ind02 = ind02(1)+1;
Phase1_offset = str2num(p(ind02:ind0));

ind0 = strfind(f2,'##$PVM_SPackArrPhase2Offset');
p = f2(ind0+4:end);
ind0 = strfind(p,')');
i1 = ind0(1)+1;
ind0 = strfind(p,'##$');
ind0 = ind0(1)-1;
ind02 = strfind(p,')');
ind02 = ind02(1)+1;
Phase2_offset = str2num(p(ind02:ind0));

ind0 = strfind(f2,'##$PVM_SPackArrSliceOffset');
p = f2(ind0+4:end);
ind0 = strfind(p,'##$');
ind0 = ind0(1)-1;
ind02 = strfind(p,')');
ind02 = ind02(1)+1;
Slice_offset = str2num(p(ind02:ind0));

ind1 = strfind(f2,'##$EffectiveTE');
if (isempty(ind1))
    %single echo
    ind1 = strfind(f2,'##$PVM_EchoTime');
    if (isempty(ind1))
        fprintf('sorry,Can not find "%s" Value,Please check Method file.\n','PVM_EchoTime');
        return;
    end
    NumEcho =1;
    p = f2(ind1+4:end);
    ind1 = strfind(p,'=');
    i1 = ind1(1)+1;
    ind1 = strfind(p,'##$');
    i2 = ind1(1)-1;
    TE = str2num(p(i1:i2))*1e-3;%should change for multiecho.
    delta_TE = TE;
else
    %multi echo
    p = f2(ind1+4:end);
    ind1 = strfind(p,'(');
    i1 = ind1(1)+1;
    ind1 = strfind(p,')');
    i2 = ind1(1)-1;
    NumEcho = str2num(p(i1:i2));
    
    p = p(i2+2:end);
    i1 =1;
    ind1 = strfind(p,'##$');
    i2 = ind1(1)-1;
    
    TE = str2num(p(i1:i2))*1e-3;%should change for multiecho.
    delta_TE = TE(2)-TE(1);
end
%%

ind1 = strfind(f2,'##$PVM_RepetitionTime');
if (isempty(ind1))
    fprintf('sorry,Can not find "%s" Value,Please check Method file.\n','PVM_RepetitionTime');
    return;
end
p = f2(ind1+4:end);
ind1 = strfind(p,'=');
i1 = ind1(1)+1;
ind1 = strfind(p,'##$');
i2 = ind1(1)-1;
TR = str2num(p(i1:i2))*1e-3;

ind1 = strfind(f2,'##$PVM_SPackArrGradOrient');
if (isempty(ind1))
    fprintf('sorry,Can not find "%s" Value,Please check Method file.\n','PVM_SPackArrGradOrient');
    return;
end
p = f2(ind1+4:end);
ind1 = strfind(p,')');
i1 = ind1(1)+1;
ind1 = strfind(p,'##$');
i2 = ind1(1)-1;
str1 = p(i1:i2);
%elmiate the return character.
t = 1;
for k = 1:length(str1)
    if (int32(str1(k))>30)
     str(t) = str1(k);
     t=t+1;
    end
end
PVM_SPackArrGradOrient = str2num(str);
af = reshape(PVM_SPackArrGradOrient,3,3);
%af(3,:) = -af(3,:);
%Affine3D = [af(:,2),af(:,1),af(:,3)];
Affine3D = af;
B0_dir = Affine3D\[0 0 1]';

ind1 = strfind(f1,'##$BF8');
if (isempty(ind1))
    fprintf('sorry,Can not find "%s" Value,Please check Method file.\n','BF8');
    return;
end
p = f1(ind1+4:end);
ind1 = strfind(p,'=');
i1 = ind1(1)+1;
ind1 = strfind(p,'##$');
i2 = ind1(1)-1;
CF = str2num(p(i1:i2))*1e6;


save parameter voxel_size matrix_size TE delta_TE CF Affine3D B0_dir TR NumEcho;
%%
xsize = matrix_size(1);
ysize = matrix_size(2);
zsize = matrix_size(3);
echo_no = NumEcho;
iField = zeros(xsize,ysize,zsize,echo_no);
%for k = 1:length(acq_no)
    raw_fid_filename=[fid_dir '/fid'];
    fid = fopen(raw_fid_filename,'r');
    data = fread(fid,Inf,'int32');
     fclose(fid);
    if xsize<=128
        xres = 128;
    elseif xsize<=256
        xres = 256;
    elseif xsize<=384
        xres = 384;
    elseif xsize<=512
        xres = 512;
    else
        xres = 1024;
    end
    data=reshape(data,[2 xres echo_no ysize zsize ]);
    data=data(:,1:xsize,:,:,:);
    data0(:,:,:,:) = complex(data(1,:,:,:,:),data(2,:,:,:,:));
   
    clear('data');
    if echo_no >1
        for k=1 : echo_no
            data1 = fftshift(ifftn(fftshift(squeeze(data0(:,k,:,:)))));
            iField(:,:,:,k) = data1;
        end
    else
        data1 = fftshift(ifftn(fftshift(squeeze(data0(:,:,:)))));
        iField(:,:,:) = data1;
    end
    
    ksp_shift = [Phase2_offset, Phase1_offset, -Slice_offset];
    FOV = matrix_size.*voxel_size;
    
    phasor = zeros(matrix_size);
    for i = 1 : matrix_size(1)
        for j = 1 : matrix_size(2)
            for k = 1 : matrix_size(3)
                phasor(i,j,k) = exp(-1i*2*pi*((i-matrix_size(1)/2)*ksp_shift(1)/FOV(1)+(j-matrix_size(2)/2)*ksp_shift(2)/FOV(2) + (k - matrix_size(3)/2)*ksp_shift(3)/FOV(3)));
            end
        end
    end
    
    for i = 1 : echo_no
        iField(:,:,:,i) = ifftn(ifftshift(fftshift(fftn(squeeze(iField(:,:,:,i)))).*phasor));
    end
    iField = conj(iField);%fit GE data Mode
    
    
    
end