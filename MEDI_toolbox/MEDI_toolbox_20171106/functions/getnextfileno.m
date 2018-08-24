% Generate file numbers for saving results
%   fileno=getnextfileno(folder, prefix, suffix)
% 
%   output
%   fileno - the next available file number
%
%   input
%   folder - the target directory
%   prefix - the prefix of the filename
%   suffix - the suffix of the filename
%
%   Created by Ildar Khalidov in 2010
%   Last modified by Tian Liu on 2013.07.24

function fileno=getnextfileno(folder, prefix, suffix)

filenames=dir(folder);
N_files=length(filenames);

nmax=0;
for k=1:N_files
    tokens=regexp(filenames(k).name,strcat('^',regexptranslate('escape',prefix),'([0-9]+)',regexptranslate('escape',suffix),'$'),'tokens');
    if(~isempty(tokens))
        nmax=max(str2num(char(tokens{1})),nmax);
    end
end

fileno=nmax+1;