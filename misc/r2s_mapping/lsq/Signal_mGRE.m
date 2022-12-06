%% S = Signal_mGRE(m0,x2s,te,varargin)
%
% Description: Multiple exponential decay model for multiecho GRE
%
% Input
% ----------
%   m0      : T1weighted signal
%   x2s     : either T2* or R2*
%   te      : echo times
%   ***flag***
%   ----------
%   'r2s'   : input R2* instead of T2*
%   'freq'  : frequency shift of all components (optional)
%
% Output
% ----------
%   s       : T2*-weighted signal
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 23 February 2017
% Date last modified: 13 October 2017
%
function s = Signal_mGRE(s0,x2s,te,varargin)
%% check and parse input
[mode,freq,verbose] = parse_varargin_Signal_mGRE(varargin);

% number of exponential componets
ncomp = length(s0);

% check if the size of m0 and x2s match
if length(x2s)~=ncomp
    error('Mismatch of m0 size of relaxation constant size');
end
% check if frequency ok
if length(freq)~=ncomp && freq(1)~=0
    error('Mismatch of frequency size of relaxation constant size');
end
if freq(1)==0 && length(freq)==1
    freq = repmat(freq(1),1,ncomp);
end
% convert t2s to r2s
if strcmpi(mode,'r2s')
    t2s = 1./x2s;
else
    t2s = x2s;
end
if verbose
    fprintf('%i exponential component(s)\n',ncomp);
end

%% Core
s = zeros(1,length(te));
for kt=1:length(te)
    s(kt) = sum(s0(:).*exp(-te(kt)./t2s(:)).*exp(1i*2*pi*freq(:)*te(kt))) ;
end

end

%% Parsing varargin
function [mode,freq,verbose] = parse_varargin_Signal_mGRE(arg)
mode = 't2s';
freq = 0;
verbose=false;
for kvar = 1:length(arg)
    if strcmpi(arg{kvar},'r2s')
        mode = 'r2s';
        continue
    end
    if strcmpi(arg{kvar},'freq')
        freq = arg{kvar+1};
        continue
    end
    if strcmpi(arg{kvar},'-v')
        verbose = true;
        continue
    end
end
end