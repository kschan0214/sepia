function  runmex(  )
%RUNMEX Summary of this function goes here
%   Detailed explanation goes here
opt = {['-I.'],...
       ['-I' fullfile('fsl','extras','include','newmat')], ...
       ['-I' fullfile('fsl','include','meshclass')], ...
       ['-I' fullfile('fsl','extras','include','boost')], ...
       ['-I' fullfile('fsl','include','utils')], ...
       ['-I' fullfile('fsl','include','newimage')], ...
       ['-I' fullfile('fsl','include')], ...
       ['-DEXPOSE_TREACHEROUS']};

srcf = 'bet2.cpp';
tmp = dir('*.c*');
sources = {srcf};
n=0;
for i=1:length(tmp);
    if ~strcmp(tmp(i).name,srcf)
        n = n + 1;
        sources{n+1}=tmp(i).name;
    end
end
mex('-v',opt{:},sources{:})

end

