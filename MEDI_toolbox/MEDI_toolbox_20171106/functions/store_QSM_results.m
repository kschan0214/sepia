% An interal function to parse input arguments for MEDI_L1
%   Save the results in a certain folder
%   Default folder is ./results/
%   Adapted from Ildar Khalidov
%   Modified by Tian Liu on 2011.02.02
%   Last modified by Tian Liu on 2013.07.24


function store_QSM_results(varargin )
[QSM, summary, iMag, RDF, Mask] = parse_input(varargin{:});

if exist('results','dir') == 0
    mkdir results;
end

fileno=getnextfileno('results/','x','.mat');
save(strcat('results/x',sprintf('%08u',fileno), '.mat'), 'QSM', 'summary','iMag','RDF','Mask');
end

function [QSM, summary, iMag, RDF, Mask] = parse_input(varargin)
QSM = varargin{1};
iMag = varargin{2};
RDF = varargin{3};
Mask = varargin{4};
summary = cell2struct(varargin(6:2:end),varargin(5:2:end),2 );
end
