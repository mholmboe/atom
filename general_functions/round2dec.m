%% round2dec.m
% * This function rounds a value to a certain decimal, and assures
% compatibility between MATLAB and Octave which can't do round2dec(X,precision)
%
%% Version
% 3.00
%
%% Contact
% Please report problems/bugs to michael.holmboe@umu.se
%
%% Examples
% # out_value = round2dec(in_value,varargin)
%

function out_value = round2dec(in_value,varargin)

if nargin==1
    out_value = round(in_value);
elseif nargin>1
    precision=varargin{1};
    out_value = round2dec(in_value.*(10^precision))./(10^precision);
else
    disp('Could not find input variable')
end