% Welcome to STI SUITE!
%
% The STI Suite software package was developed by:
% Wei Li, PhD and Chunlei Liu, PhD
% Release 1.0   08/29/2013
%         1.1   12/01/2013, minor correction to V_SHARP_v1.m
%               12/11/2013, corrected image loading buttons "mat" for
%               magnitude, phase and mask.
%         1.4   01/22/2014, Major update. Corrected code for GUIs, allow 
%               zooming, and removed errnonous depedencies;
%               User Manual updated.
%         2.0   03/17/2014, Major update. Wei Li
%               Better Brain Extraction, HARPERELLA without change of mask 
%               and others;
%               User Manual updated.
%         2.01   6/15/2014, Major update. Wei Li
%               Better methods for HARPERELLA, QSM, Upgraded QSM_GUI.
%               Working version name: STI Vista.
%         2.2   01/05/2015, Minor update. Wei Li
%               Odd number of voxels.
%               Removed the expiration date for the tempratory version.
%               All data were converted into single format. If you need
%               double precision, please feel free to contact us.
% Chunlei Liu, Duke University.
% Wei Li, Duke University, University of Texas Health Science Center.
% (c) All rights reserved.

% ------------------------------TERMS OF USE-------------------------------
%  
% The algrithms, including the Laplacian-based processing and quantative 
% susceptibility mapping methods, are free for academic use. Any commerical
% or industrial use is not permitted without the permission of the authors.
%
% We will appreciate it if you could cite our relevant publications.

% Laplacian-based phase processing (phase unwrapping, V-SHARP, HARPERELLA)
% [1] Li W, Wu B, Liu C. Quantitative susceptibility mapping of human 
% brain reflects spatial variation in tissue composition, NeuroImage. 2011; 
% 15; 55:1645.
% [2] Li W, Wu B, Avram AV, Liu C. Integrated Phase Unwrapping and 
% Background Phase Removal using Laplacian for Quantitative Susceptibility 
% Mapping. NMR Biomed 2014,27: 219-227
% [3] Wu B, Li W, Liu C. Whole brain susceptibility mapping using 
% compressed sensing. Mag. Res. Med. 2012;67(1):137. (BW and WL contributed
% equally).

% The LSQR and iLSQR methods for quantitative susceptibility mapping (QSM)
% [4] Li W, Wu B, Liu C. Quantitative susceptibility mapping of human 
% brain reflects spatial variation in tissue composition, NeuroImage. 2011; 
% 15; 55:1645.
% [5] Li W, Wang N, Yu F, Han H, Cao W, Romero R, Tantiwongkosi B, Duong TQ, 
% Liu C. A method for estimating and removing streaking artifacts in 
% quantitative susceptibility mapping. NeuroImage, 2015;108:111–122

% The Susceptibility Tensor Imaging (STI)
% [6] Liu C. Susceptibility Tensor Imaging. Magn Reson Med. 2010;63(6):
% 1471-7
% [7] Li W, Wu B, Avram AV, Liu C, Magnetic susceptibility anisotropy of in
% vivo human brain and its molecular underpinnings. NeuroImage. 2012;59(3):
% 2088.
% [8] Liu C, Li W, Wu B, Jiang Y, Johnson GA. 3D fiber tractography with 
% susceptibility tensor imaging. NeuroImage. 2012;59(2):1290.
% [9] Li W, Liu C. Comparison of Magnetic Susceptibility Tensor and 
% Diffusion Tensor of the Brain. Journal of Neuroscience and 
% Neuroengineering. 2014;2;431-440
%
% The GUI tools is designed for free academic use, please feel free to 
% modify it for your applications.
% 
% For any commerical or industrial usage, please contact Dr Chunlei Liu at
% Chunlei.Liu@duke.edu
% 
% We provide support at our discretion (no guarranty, sorry).
% For technical questions, please contact STI.Suite.MRI@gmail.com
%
%



