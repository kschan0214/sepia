% An interal function to parse input arguments for MEDI_L1
%   Created by Tian Liu on 2011.03.16
%   Modified by Tian Liu and Shuai Wang on 2011.03.15
%   Last modified by Tian Liu on 2013.07.24
% Date last modified: 27 Feb 2020 (v0.8.0, KC)

function [N_std,iMag,lam,pad,delta_TE,CF,B0_dir,merit,smv,radius,data_weighting,gradient_weighting,Debug_Mode,lam_CSF,Mask_CSF,tmp_output_dir,percentage] = parse_varargin_MEDI_L1(arg)
gyro = 42.57747892;

B0_dir              = [0,0,1];
CF                  = gyro*1e6*3;    %3T
merit               = 0;
smv                 = 0;
lam                 = 1000;
radius              = 5;
data_weighting      = 1;
gradient_weighting  = 1;
pad                 = 0;
Debug_Mode          = 'NoDebug';
percentage          = 90;

% CSF regularization
lam_CSF = 100;
Mask_CSF = [];
N_std = [];
iMag = [];
tmp_output_dir = [pwd filesep];

if ~isempty(arg)
    for kvar = 1:length(arg)
        if strcmpi(arg{kvar},'noisestd')
            N_std=arg{kvar+1};
        end
        if strcmpi(arg{kvar},'magnitude')
            iMag=arg{kvar+1};
            if size(iMag,4) > 1
                iMag = sqrt(sum(abs(iMag).^2,4));
            end
        end
        if strcmpi(arg{kvar},'lambda')
            lam=arg{kvar+1};
        end
        if strcmpi(arg{kvar},'data_weighting')
            data_weighting=arg{kvar+1};
        end
        if strcmpi(arg{kvar},'gradient_weighting')
            gradient_weighting=arg{kvar+1};
        end
        if strcmpi(arg{kvar},'merit')
            merit=arg{kvar+1};
        end
        if strcmpi(arg{kvar},'smv')
            smv = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'radius')
            radius = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'zeropad')
            pad = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'DEBUG')
            Debug_Mode = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'lambda_CSF')
            lam_CSF = arg{kvar+1};
        end        
        % CSF regularization
        if strcmpi(arg{kvar},'Mask_CSF')
            Mask_CSF = arg{kvar+1};
        end 
        if strcmpi(arg{kvar},'TE')
            delta_TE = arg{kvar+1};
        end 
        if strcmpi(arg{kvar},'b0dir')
            B0_dir = arg{kvar+1};
        end 
        if strcmpi(arg{kvar},'CF')
            CF = arg{kvar+1};
        end 
        if strcmpi(arg{kvar},'tmp_output_dir')
            tmp_output_dir = arg{kvar+1};
        end
        if strcmpi(arg{kvar},'percentage')
            percentage = arg{kvar+1};
        end
    end
end
end
