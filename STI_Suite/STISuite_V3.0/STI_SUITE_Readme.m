% Welcome to STI SUITE!
%
% The STI Suite software package was developed by:
% Hongjiang Wei, PhD, Wei Li, PhD, and Chunlei Liu, PhD
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
%         3.0   1/05/2017, Major update. Hongjiang Wei
%               Better methods STAR for QSM, Upgraded QSM_GUI.
%               2D phase processing for QSM based on 2D EPI data 
%               Working version name: STISuite_V3.0.
% Chunlei Liu, University of California, Berkeley.
%
% (c) All rights reserved.

% ------------------------------TERMS OF USE-------------------------------
%  
% The algrithms, including the Laplacian-based processing and quantative 
% susceptibility mapping methods, are free for academic use. Any commerical
% or industrial use is not permitted without the permission of the authors.
%
% We will appreciate it if you could cite our relevant publications.

%% STAR-QSM
% [1] Wei H, Dibb R, Zhou Y, Sun Y, Xu J, Wang N, Liu C. Streaking artifact
% reduction for quantitative susceptibility mapping of sources with large 
% dynamic range. NMR Biomed, 2015;28;1294-303
% [2] Wei H, Xie L, Dibb R, Li W, Decker K, Zhang Y, Johnson A, Liu C. 
% Streaking artifact reduction for quantitative susceptibility mapping of 
% sources with large dynamic range. NMR Biomed, 2015;28;1294-303
% [3] Wei H, Dibb R,  Decker K, Wang Nian, Zhang Y, Zong X, Nissman D, Liu C.
% Investigating Magnetic Susceptibility of Human Knee Joint at 7 Tesla. 
% MRM [Epub ahead of print]

%% 2D phase Processing for QSM
% [4] Wei H, Zhang Y, Gibbs E,Chen N,Wang N, Liu C. Joint 2D and 3D phase 
% processing for quantitative susceptibility mapping: application to 2D 
% echo-planar imaging. 2016 Feb 17. doi: 10.1002/nbm.3501. [Epub ahead of print]

%% Laplacian-based phase processing (phase unwrapping, V-SHARP)
% [5] Li W, Wu B, Liu C. Quantitative susceptibility mapping of human 
% brain reflects spatial variation in tissue composition, NeuroImage. 2011; 
% 15; 55:1645.
% [6] Wu B, Li W, Liu C. Whole brain susceptibility mapping using 
% compressed sensing. Mag. Res. Med. 2012;67(1):137. (BW and WL contributed
% equally).

% The LSQR method for quantitative susceptibility mapping (QSM)
% [7] Li W, Wu B, Liu C. Quantitative susceptibility mapping of human 
% brain reflects spatial variation in tissue composition, NeuroImage. 2011; 
% 15; 55:1645.

% The Susceptibility Tensor Imaging (STI)
% [8] Liu C. Susceptibility Tensor Imaging. Magn Reson Med. 2010;63(6):
% 1471-7
% [9] Li W, Wu B, Avram AV, Liu C, Magnetic susceptibility anisotropy of in
% vivo human brain and its molecular underpinnings. NeuroImage. 2012;59(3):
% 2088.
% [10] Liu C, Li W, Wu B, Jiang Y, Johnson GA. 3D fiber tractography with 
% susceptibility tensor imaging. NeuroImage. 2012;59(2):1290.


% The GUI tools is designed for free academic use, please feel free to 
% modify it for your applications.
% 
% For any commerical or industrial usage, please contact Dr Chunlei Liu at
% chunlei.liu@berkeley.edu
% 
% We provide support at our discretion (no guarranty, sorry).
% For technical questions, please contact STI.Suite.MRI@gmail.com
%
%



