function [nii] = load_nii_xfm(filename, img_idx, old_RGB, preferredForm)

   if ~exist('filename','var'),
      error('Usage: [nii] = load_nii(filename, [img_idx], [old_RGB])');
   end
   
   if ~exist('img_idx','var'), img_idx = []; end
   if ~exist('old_RGB','var'), old_RGB = 0; end
   if ~exist('preferredForm','var'), preferredForm= 's'; end     % Jeff
   
   %  Read the dataset header
   %
   [nii.hdr,nii.filetype,nii.fileprefix,nii.machine] = load_nii_hdr(filename);

   %  Read the header extension
   %
%   nii.ext = load_nii_ext(filename);
   
   %  Read the dataset body
   %
   [nii.img,nii.hdr] = ...
        load_nii_img(nii.hdr,nii.filetype,nii.fileprefix,nii.machine,img_idx,'','','',old_RGB);
   
   %  Perform some of sform/qform transform
   %
%   nii = xform_nii(nii, preferredForm);


   hdr = nii.hdr;

   %  NIFTI can have both sform and qform transform. This program
   %  will check sform_code prior to qform_code by default.
   %
   %  If user specifys "preferredForm", user can then choose the
   %  priority.					- Jeff
   %
   useForm=[];					% Jeff

   if isequal(preferredForm,'S')
       if isequal(hdr.hist.sform_code,0)
           error('User requires sform, sform not set in header');
       else
           useForm='s';
       end
   end						% Jeff

   if isequal(preferredForm,'Q')
       if isequal(hdr.hist.qform_code,0)
           error('User requires sform, sform not set in header');
       else
           useForm='q';
       end
   end						% Jeff

   if isequal(preferredForm,'s')
       if hdr.hist.sform_code > 0
           useForm='s';
       elseif hdr.hist.qform_code > 0
           useForm='q';
       end
   end						% Jeff
   
   if isequal(preferredForm,'q')
       if hdr.hist.qform_code > 0
           useForm='q';
       elseif hdr.hist.sform_code > 0
           useForm='s';
       end
   end						% Jeff

   if isequal(useForm,'s')
      R = [hdr.hist.srow_x(1:3)
           hdr.hist.srow_y(1:3)
           hdr.hist.srow_z(1:3)];

      T = [hdr.hist.srow_x(4)
           hdr.hist.srow_y(4)
           hdr.hist.srow_z(4)];

      nii.hdr.hist.old_affine = [ [R;[0 0 0]] [T;1] ];

   elseif isequal(useForm,'q')
      b = hdr.hist.quatern_b;
      c = hdr.hist.quatern_c;
      d = hdr.hist.quatern_d;

      if 1.0-(b*b+c*c+d*d) < 0
         if abs(1.0-(b*b+c*c+d*d)) < 1e-5
            a = 0;
         else
            error('Incorrect quaternion values in this NIFTI data.');
         end
      else
         a = sqrt(1.0-(b*b+c*c+d*d));
      end

      qfac = hdr.dime.pixdim(1);
      i = hdr.dime.pixdim(2);
      j = hdr.dime.pixdim(3);
      k = qfac * hdr.dime.pixdim(4);

      R = [a*a+b*b-c*c-d*d     2*b*c-2*a*d        2*b*d+2*a*c
           2*b*c+2*a*d         a*a+c*c-b*b-d*d    2*c*d-2*a*b
           2*b*d-2*a*c         2*c*d+2*a*b        a*a+d*d-c*c-b*b];

      T = [hdr.hist.qoffset_x
           hdr.hist.qoffset_y
           hdr.hist.qoffset_z];

      nii.hdr.hist.old_affine = [ [R * diag([i j k]);[0 0 0]] [T;1] ];

   elseif nii.filetype == 0 & exist([nii.fileprefix '.mat'],'file')
      load([nii.fileprefix '.mat']);	% old SPM affine matrix
      R=M(1:3,1:3);
      T=M(1:3,4);
      T=R*ones(3,1)+T;
      M(1:3,4)=T;
      nii.hdr.hist.old_affine = M;

   else
      M = diag(hdr.dime.pixdim(2:5));
      M(1:3,4) = -M(1:3,1:3)*(hdr.hist.originator(1:3)-1)';
      M(4,4) = 1;
      nii.hdr.hist.old_affine = M;
   end

   return					% load_nii_no_xform
