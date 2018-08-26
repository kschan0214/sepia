function ni = readFileNifti(fileName, volumesToLoad)
%
% niftiImage = readFileNifti(fileName, [volumesToLoad=-1])
%
% Reads a NIFTI image and populates a structure that should resemble
% the NIFTI 1 standard (see http://nifti.nimh.nih.gov/nifti-1/ ).
% The optional second arg specifies which volumes to load for a 4-d dataset.
% The default (-1) means to read all, [] (empty) will just return the header.
%
% Call this function again with an output argument to get an empty structure.
%
% NOTE: this file contains a very slow m-file implementation of a compiled
% mex function. If you get a warning that the mex function is not being
% called, then compiling readFileNifti.c will dramatically improve
% performance.
%
% HISTORY:
% 2008.10.01 RFD: wrote it.
% 2009.07.01 RFD: fixed c-to-matlab 0-to-1 offset error. The qto and sto
% offsets prior to this were shifted by 1 voxel.

% disp('NOTE: Compile readFileNifti.c for faster operation!');

if(nargin==0)
    if(nargout==0)
        help(mfilename);
    else
        ni = getNiftiStruct;
        return;
    end
end
    
if(~exist('volumesToLoad','var') || volumesToLoad==-1)
    volumesToLoad = []; 
end
ni = getNiftiStruct;
ni.fname = fileName;
[p,f,e] = fileparts(fileName);
if(strcmpi(e,'.gz'))
   tmpDir = tempname;
   mkdir(tmpDir);
   gunzip(fileName, tmpDir);
   tmpFileName = fullfile(tmpDir, f);
   tmpFile = true;
   nii = load_nii(tmpFileName, volumesToLoad);
else
   tmpFile = false;
   nii = load_nii(fileName, volumesToLoad);
end
ni.data = nii.img;
ni.ndim = ndims(ni.data);
ni.dim = size(ni.data);
ni.pixdim = nii.hdr.dime.pixdim(2:ni.ndim+1);
ni.scl_slope = nii.hdr.dime.scl_slope;
ni.scl_inter = nii.hdr.dime.scl_inter;
ni.cal_min = nii.hdr.dime.cal_min;
ni.cal_max = nii.hdr.dime.cal_max;
ni.slice_code = nii.hdr.dime.slice_code;
ni.slice_start = nii.hdr.dime.slice_start;
ni.slice_end = nii.hdr.dime.slice_end;
ni.slice_duration = nii.hdr.dime.slice_duration;
ni.toffset = nii.hdr.dime.toffset;
% FIX ME: decode the units
% 'unknown,meter,mm,micron,sec,msec,usec,hz,ppm,rads,unknown'
% 'unknown,meter,mm,micron,sec,msec,usec,hz,ppm,rads,unknown'
ni.xyz_units = nii.hdr.dime.xyzt_units;
ni.time_units = nii.hdr.dime.xyzt_units;
ni.intent_code = nii.hdr.dime.intent_code;
ni.intent_p1 = nii.hdr.dime.intent_p1;
ni.intent_p2 = nii.hdr.dime.intent_p2;
ni.intent_p3 = nii.hdr.dime.intent_p3;
ni.intent_name = nii.hdr.hist.intent_name;
% CHECK THIS:
ni.nifti_type = nii.hdr.dime.datatype;
ni.qform_code = nii.hdr.hist.qform_code;
ni.sform_code = nii.hdr.hist.sform_code;
% The fields freq_dim, phase_dim, slice_dim are all squished into the single
% byte field dim_info (2 bits each, since the values for each field are
% limited to the range 0..3):
ni.freq_dim =  bitand(uint8(3),          uint8(nii.hdr.hk.dim_info)    );
ni.phase_dim = bitand(uint8(3), bitshift(uint8(nii.hdr.hk.dim_info),-2));
ni.slice_dim = bitand(uint8(3), bitshift(uint8(nii.hdr.hk.dim_info),-4));
ni.quatern_b = nii.hdr.hist.quatern_b;
ni.quatern_c = nii.hdr.hist.quatern_c;
ni.quatern_d = nii.hdr.hist.quatern_d;
ni.qoffset_x = nii.hdr.hist.qoffset_x;
ni.qoffset_y = nii.hdr.hist.qoffset_y;
ni.qoffset_z = nii.hdr.hist.qoffset_z;
ni.qfac = nii.hdr.dime.pixdim(1);
% We need to apply the 0-to-1 correction to qto and sto matrices.
offset = [ni.qoffset_x, ni.qoffset_y, ni.qoffset_z] - ni.pixdim(1:3);
ni.qto_xyz = quatToMat(ni.quatern_b, ni.quatern_c, ni.quatern_d, offset(1), offset(2), offset(3), ni.qfac, ni.pixdim);
if(all(ni.qto_xyz(1:9)==0)), ni.qto_ijk = zeros(4); 
else                         ni.qto_ijk = inv(ni.qto_xyz); end
ni.sto_xyz = [nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; [0 0 0 1]];
ni.sto_xyz(1:3,4) = ni.sto_xyz(1:3,4) - ni.pixdim(1:3)';
if(all(ni.sto_xyz(1:9)==0)), ni.sto_ijk = zeros(4); 
else                         ni.sto_ijk = inv(ni.sto_xyz); end

ni.descrip = nii.hdr.hist.descrip;
ni.aux_file = nii.hdr.hist.aux_file;

if(tmpFile)
    delete(tmpFileName);
end
return;


function ni = getNiftiStruct()
ni.data = [];
ni.fname = '';
ni.ndim = [];
ni.dim = [];
ni.pixdim = [];
ni.scl_slope = [];
ni.scl_inter = [];
ni.cal_min = [];
ni.cal_max = [];
ni.qform_code = [];
ni.sform_code = [];
ni.freq_dim = [];
ni.phase_dim = [];
ni.slice_dim = [];
ni.slice_code = [];
ni.slice_start = [];
ni.slice_end = [];
ni.slice_duration = [];
ni.quatern_b = [];
ni.quatern_c = [];
ni.quatern_d = [];
ni.qoffset_x = [];
ni.qoffset_y = [];
ni.qoffset_z = [];
ni.qfac = [];
ni.qto_xyz = [];
ni.qto_ijk = [];
ni.sto_xyz = [];
ni.sto_ijk = [];
ni.toffset = [];
ni.xyz_units = 'unknown,meter,mm,micron,sec,msec,usec,hz,ppm,rads,unknown';
ni.time_units = 'unknown,meter,mm,micron,sec,msec,usec,hz,ppm,rads,unknown';
ni.nifti_type = [];
ni.intent_code = [];
ni.intent_p1 = [];
ni.intent_p2 = [];
ni.intent_p3 = [];
ni.intent_name = [];
ni.descrip = [];
ni.aux_file = [];
return;
