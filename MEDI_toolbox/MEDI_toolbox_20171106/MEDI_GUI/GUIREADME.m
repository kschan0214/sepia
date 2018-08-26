% This is the readme file for the GUI beta 03 with stripBD version
% This readme gives the guide for how to use this gui and toolbox
% Author: Yuecong "Arch" Wu
% Email: archwu0817@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To use this GUI, add the folder of location of gui to your MATLAB Path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Loading is the first step to properly use this gui. click Browse 
% to select the directory data located at and click Select Folder button
% after navigated to the desired folder. Then click Load button to start
% loading the data.
% 
% Data Loading is needed before performing any other actions for this gui
% The data directory is required to be filled into the editable text box
% before doing the 'Save DICOM'
% This GUI supports the DICOM files from Siemens, GE and Philips.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Applying Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are two options for applying mask to the graph. 'Auto' and 
% 'User Select'. 'Auto' mode, literally automatically generates a mask for
% use. 'User Select' enables the user to select the mask as user demands.
% When Using 'User Select' Mode, the user can use the browse button to 
% select the mask file. After selecting the mask the user can click 'Mask' 
% button to load the mask into the workspace.
% 
% This GUI supports the .img and .mat file as the Mask.
% Auto mode would generate the mask right after changing selection to
% 'Auto', this action requires loading the data before
% StripBD would perform a strip for current loaded mask.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Applying Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Run Modes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% There are two running modes 'Default' and 'Advanced'. The user can 
% change running mode by changing the selection of the radio buttons.
% 'Advanced' mode will enable the modification of Unwrapping mode, 
% BF Removing mode and the input variable for MEDI.
% 
% 'Default' mode will enable the 'JUST DO IT' button which would auto-
% matically perform all the functions be default parameters once clicked.
%
% Just Do It won't perform save DICOM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Run Modes   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Functions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The user should perform the functions in correct order, but there is
% a feature that as long as the user loaded the data and mask(generated
% mask is also good) and have made selections of Unwrapping mode, 
% BFRemoval mode and parameters of MEDI. The GUI would automatically
% execute. First Fitting, then Unwrapping and BF Removing then finally
% save the RDF into a folder as '*_result' and MEDI. '*' represents the
% folder of the DICOM user have loaded in.
% 
% 'Save DICOM' Retrieves infomation from original DICOMs and QSM, then  
% save to new DICOM files into the folder'*_result/DICOM/'. '*' represents 
% the folder of the DICOM user have loaded in.
%
% The default value for MEDI is lambda - 1000 edge - 0.90 SMV - 5(mm)
% Leaving the Edge and SMV editable text blank or with value that is not
% a valid number would result in running MEDI with default value as above
%
% Please refer to README.m for detailed information of functions used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Functions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%    Multi-Process   %%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-datasets mode enables the user to process multiple dataset at a
% click after manually added all the datapaths the user desires to
% process. By clicking 'MultiBrowse' button the user can add data 
% directories to the listbox.
% If the user want to clear the selection, 'Clear Directories' may do the
% job.
% After making selections the user may click on 'MultiProcess' and the GUI
% will start processing and it may take a long time depends on the size
% of datasets.
% 'MultiProcess' performs the processing using auto generates mask and the 
% user defined parameters in the left panels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%    Multi-Process   %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%   Variable inspector   %%%%%%%%%%%%%%%%%%%%%%%%%
% The panel in the right is a part for users to verify the data and
% results of calculations. Press 'Update Listbox' can show all the values
% in workspace and the user can make selection and click 'Vis' button to
% visualize the data to check the validity.
% This is supposed to be used under single dataset step-by-step using
%%%%%%%%%%%%%%%%%%%%%%%%%   Variable inspector   %%%%%%%%%%%%%%%%%%%%%%%