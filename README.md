# qsm_hub (work in progress)
==================================================

## Introduction
--------------------------------------------------

Welcome to the beta version of **QSM Hub**. **QSM Hub** is a tool providing a graphical user 
interface to build data processing pipeline of quantitative susceptibility mapping (QSM) in Matlab.

This GUI is built based on three toolboxes including [MEDI](http://weill.cornell.edu/mri/pages/qsm.html) 
(version 2017-11-06), [STI Suite](https://people.eecs.berkeley.edu/~chunlei.liu/software.html) (version 3.0)
and [FANSI](https://gitlab.com/cmilovic/FANSI-toolbox).

**QSM Hub** provides two key features for QSM processing:

1. mix-and-match methods from different toolboxes to build your own QSM processing pipeline
2. graphical user interface to easily adjust parameters of different algorithms

The objective of **QSM Hub** is to provide a platform for easy access of different QSM processing 
method in the field. To achieve this, most of the codes were written for data flow control and 
algorithm parameter control. Through **QSM Hub** I hope researchers which are not coming from the field 
would also be able to use QSM for their research. For those more experienced researchers familiar with 
QSM, I also hope this tool can simplify your workflow. 

This toolbox is still under development, so you may encounter some bugs.

**QSM Hub** is freely available for educational and research purposes. Commerical or clinical use is not 
permitted. When you use **QSM Hub** in your research, please cite the related papers of the methods 
in your processing pipeline. Details of each processing method is provided with the hyperlinks in this
document. Clicking the hyperlinks of the methods below will redirect you to the 
paper webiste or you can also find the full details in the Referencing section below. 

If you have any question or you would like to provide suggestion to improve this toolbox/report 
bug(s) please feel free to contact me k.chan@donders.ru.nl (Kwok-Shing Chan).

Kwok  
2018-05-30

----------------------------------------------------------------------------------------------------

## Installation
--------------------------------------------------

To use **QSM Hub**, add the directory containing qsm_hub.m into your MATLAB PATH

This can be done by:  
'Set Path' -> 'Add Folder' -> /your/qsm_hub/directory/ -> 'Save'  
  
or
  
with MATLAB's command: `addpath('/your/qsm_hub/directory')`  

Then you have to specify the directories of each toolbox in SpecifyToolboxesDirectory.m:

``  
MEDI_dir = '/path/to/MEDI/toolbox/';
STISuite_dir = '/path/to/STISuite/toolbox/';
FANSI_dir = '/path/to/FASI/toolbox/';  
``

Then, you can start the GUI by entering `qsm_hub` in the MATLAB's command window.  

----------------------------------------------------------------------------------------------------
