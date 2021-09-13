function print_par(par,varargin)

% This little function prints a parameter structure or matrix into a form
% that can be directly copied into your script. This can be handy, for
% example when the profiling reports a better optimum.
%
% Author     : Tjalling Jager
% Date       : February 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo2

names = glo2.names; % extract names structure from global

if numel(varargin)>0
    fid = varargin{1}; % file to write to OR flag for displaying fitted parameters
else
    fid = -1; % otherwise, print to screen
end

if isstruct(par) % if input is a structure, first unpack it
    pmat = packunpack(1,par,0);  % transform structure into a regular matrix
else
    pmat = par; % otherwise, it already is a regular matrix
end
nfields = length(names); % how many parameters do we have

% if fid == -1 % print on screen
%     disp(' ')
%     disp('Lines below may be directly copied into your script')
%     disp(' ')
% end

% Now print the matrix in a format that we want
for i = 1:nfields
    switch fid
        case -1 % print all parameters on screen
            fprintf('par.%-6s = [%10.5g  %1.0f %10.5g %10.5g %1.0f]; \n',names{i},pmat(i,1:5));
        case -2 % only print fitted parameters
            if pmat(i,2) == 1
                fprintf('par.%-6s = [%10.5g  %1.0f %10.5g %10.5g %1.0f]; \n',names{i},pmat(i,1:5));
            end
        otherwise % print to a named file
            fprintf(fid,'par.%-6s = [%10.5g  %1.0f %10.5g %10.5g %1.0f]; \n',names{i},pmat(i,1:5));
    end
end
% if fid == -1 % print on screen
%     disp(' ')
% end
