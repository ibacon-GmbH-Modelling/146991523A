function [data_all,scen_tot,unit,error_flag] = load_data_openguts(filename,filepath)

% Function to load openGUTS input files for calibration, and translate them
% into the BYOM format. The main challenge was to get this working in BYOM
% when entering multiple files to analyse simultaneously. openGUTS treats
% all data sets as isolated, so treatments can have the same identifier in
% two data sets without clashing. However, in BYOM, there is one <X0mat> to
% run all scenarios, so the identifiers must be unique (if the treatment is
% different).
% 
% Author     : Tjalling Jager 
% Date       : March 2020
% Web support: http://www.debtox.info/byom.html

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <load_data> code that is
% distributed as part of the Matlab version of openGUTS (see
% http://www.openguts.info). Therefore, this code is distributed under
% the same license as openGUTS (GPLv3). The modifications are not in the
% algorithm itself but only to ensure that the code operates in the general
% BYOM framework.
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%  
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

global glo

if ~iscell(filename) % the GUI makes it a cell array only when multiple files are selected ...
    filename = {filename}; % but we also want a cell array if it's only one
end

data_all   = []; % initialise output as empty as errors prompt a premature return
data       = cell(length(filename),1); % initialise empty cell array
error_str1 = {}; % define error string as empty cell array
error_flag = 0; % unless something happens, error flag is zero

for i_d = 1:length(filename) % run through filenames that are entered
    
    % Note: if the Matlab GUI is used for selecting files, there is no
    % need to check if the file exists, I guess.
    if exist([filepath,filename{i_d}],'file') == 2
        % Note: the <text_parser> (sub-function below) is used to parse
        % the text and convert it to a structured cell array <data>.
        [data{i_d},error_str] = text_parser(filename{i_d},filepath);
        error_str1 = cat(1,error_str1,error_str); % add the error to the previous ones
    else % if the file is not available in the right folder, produce an error
        error_str1 = cat(1,error_str1,['There is no data set with filename ',filename{i_d},' ']);
    end
    
end

for i_d = 1:length(filename) % run through filenames that are entered
    if size(data{i_d}.scen,1) == 1 % then we have constant exposure concentration
        scen{i_d} = data{i_d}.scen(2:end);
        cons(i_d) = 1;
    else
        scen{i_d} = 1:length(data{i_d}.name);
        cons(i_d) = 0;
        ind_ctrl = strcmpi(data{i_d}.name,'control');
        scen{i_d}(ind_ctrl) = 0; % call that one identifier zero!
    end
end

% for the data sets with constant exposure, we can combine all exposure
% indentifiers and keep the unique ones
ind_cons = find(cons == 1); % which data sets have constant exposure?
scen_tot = [];
for i = 1:length(ind_cons)
    scen_tot = cat(2,scen_tot,scen{i});
end
if ~isempty(scen_tot)
    scen_tot = unique(scen_tot);
end

% For time-varying exposure, we need to make sure that the identifiers do
% not duplicate any identifier of the constant exposures, nor of other
% time-varying data sets! This requires some trickery ... also to keep the
% zeros in there, as they SHOULD be duplicated.
ind_var = find(cons == 0); % which data sets have constant exposure?
for i = 1:length(ind_var)
    scen_tmp = scen{ind_var(i)};
    ind_zero = (scen_tmp==0);
    scen_tmp(ind_zero) = []; % zeros can be duplicates
    % the other scenarios can NOT be the same as any other that we
    % currently have!
    while ~all(ismember(scen_tmp,scen_tot)==0)
        scen_tmp = scen_tmp + 10; % just add 10 to the identifiers
    end
    scen{ind_var(i)} = zeros(size(scen{ind_var(i)}));
    scen{ind_var(i)}(~ind_zero) = scen_tmp;
    scen_tot = cat(2,scen_tot,scen_tmp);
end

scen_tot = unique([0 scen_tot]); % make sure there is a control in there ...
% This is a bit silly, but the control may be removed in the previous for
% loop.

for i_d = 1:length(filename) % run through filenames that are entered

    if size(data{i_d}.scen,1) == 1 % then we have constant exposure concentration
        data_all_tmp  = [-1 data{i_d}.scen(2:end);data{i_d}.surv ]; % note that first column in scen is time (here only a zero)
        
    else
        
        scen_data = data{i_d}.scen;
%         if length(unique(scen_data(:,1))) ~= length(scen_data(:,1))
%             % if we have a renewals scenario, with time points occurring
%             % twice, extract a tiny bit from the first time point
%             scen_data(diff(scen_data(:,1))==0,1) = scen_data(diff(scen_data(:,1))==0,1) - 0.01;
%             % error_str1 = cat(1,error_str1,['For now, double time points in the exposure scenario not allowed (',filename{i_d},')']);
%         end

        % Note: adding a tiny amount to the time vector is no longer
        % needed, as make_scen is adapted to allow double time points as
        % well.
        
        % Note that scen contains the truly unique identifiers for the
        % time-varying exposure treatments (plus zero for the control)
        Cw = [1 scen{i_d};scen_data];
        Cw(:,Cw(1,:)==0) = []; % remove the controls (identifier zero) as we don't make a scenario for those
        
        data_all_tmp = [-1 scen{i_d};data{i_d}.surv]; % note that first column in scen is time
        
        % create entries for the LabelTable
        Scenario = scen{i_d}';
        Label    = data{i_d}.name';
        glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels
        
        make_scen(4,Cw); % create the globals to define the forcing function
    end
    
    if length(filename) == 1
        data_all{1} = data_all_tmp;
    else
        data_all{i_d,1} = data_all_tmp;
    end
    
end

unit = data{i_d}.concunit;

if ~isempty(error_str1)
    f = errordlg(error_str1,'Errors in input file'); % display a message box with the errors
    error_flag = 1; % signal main script that we ran into errors
    return
end

if length(filename) == 1
    glo.basenm = [glo.basenm,'_',filename{1}(1:end-4)]; % change the basename to include the FIRST filename
else
    glo.basenm = [glo.basenm,'_multi_',filename{1}(1:end-4)]; % change the basename to include the FIRST filename
end

disp('Following data sets are loaded for calibration:')
for i_d = 1:length(filename) % run through data sets loaded
    disp(['    set: ',num2str(i_d),', file: ',filename{i_d}]) % and display their filename
end

%% sub-function text_parser

function [data,error_str] = text_parser(filename,filepath)

% A sub-function to read the formatted text file for the input data. It
% reads one line at a time and interprets is. It returns the structured
% cell array <data>. 
%
% Outputs
% <data>            cell array, with one structured block per data set. 
% fields in <data>, for each data set:
%   <data.name>       names of the treatments as cell array of strings
%   <data.surv>       survival data set without headers (matrix, with time in first column)
%   <data.scen>       exposure scenario data without headers (matrix, with time in first column)
%   <data.concunit>   concentration unit of the exposure data set.
% 
% <error_str>  an array with strings of error messages

fileID    = fopen([filepath,filename]); % open file for reading
nxtline   = fgetl(fileID);   % reads one line of text
error_str = {}; % define error string as empty cell array
data      = []; % initialise as empty as errors prompt a premature return

while isempty(strfind(nxtline,'Survival time')) && isempty(strfind(nxtline,'survival time'))
% This is not the recommended syntax for Matlab. However, <contains> was
% added in release 2016b, so this formulation works error free in older
% versions. Preferred syntax for newer versions:
% % while ~contains(nxtline,'Survival time','IgnoreCase',true) % continue until we find the line that starts with "Survival time"
    nxtline = fgetl(fileID); % reads one line of text
     if nxtline == -1
         error_str = sprintf('Format error in data set %1s: unexpected end of file',filename);
         return
     end
end
surv_header = strsplit(nxtline,'\t','CollapseDelimiters',0); % split text string at the tab and makes cell array
% Multiple tabs are treated as multiple tabs (hopefully, this catches empty
% cells as well).

if length(surv_header) ~= length(unique(surv_header))
    error_str = sprintf('Format error in data set %1s: treatment names are not unique',filename);
    return
end
if sum(strcmpi(surv_header,'Control')) > 1 % is this already caught by prvious check?
    % Note: <strcmpi> compares strings without caring for case.
    error_str = sprintf('Format error in data set %1s: there can be only one treatment identified as control',filename);
    return
end

nxtline    = fgetl(fileID); % reads the next line from the file
surv_nmbrs = []; % initialise an empty matrix for survivors

while isempty(strfind(nxtline,'Concentration')) && isempty(strfind(nxtline,'concentration')) % continue until we find the line that starts with "Concentration"
% This is not the recommended syntax for Matlab. However, <contains> was
% added in release 2016b, so this formulation works error free in older
% versions. Preferred syntax for newer versions:
% % while ~contains(nxtline,'Concentration','IgnoreCase',true) % continue until we find the line that starts with "Concentration"
    
    cell_line = strsplit(nxtline,'\t','CollapseDelimiters',0); % return cell array with expressions, split at the tabs
    if length(surv_header) ~= length(cell_line)
        error_str = sprintf('Format error in data set %1s: survival data does not have same number of columns at each time point',filename);
        return
    end
    
    surv_nmbrs = [surv_nmbrs; nan(1,length(cell_line))]; % initialise a line of NaNs
    for i = 1:length(cell_line) % run through all cells of the array
        surv_i = str2double(cell_line{i}); % convert string to number
        surv_nmbrs(end,i) = surv_i; % note that any text fields will automatically become NaN
    end
    nxtline = fgetl(fileID); % reads next line from file
    if nxtline == -1
        error_str = sprintf('Format error in data set %1s: unexpected end of file',filename);
        return
    end
end 

conc_unit = strsplit(nxtline,'\t'); % split text string at the tab and makes cell array
% The initial version of the openGUTS data format applied a space between
% the words "Concentration unit:" and the unit itself. We modified it for
% the software to have a tab in between, which is more logical. However, I
% here want to be able to use the old version as well. 
if numel(conc_unit) == 1 || isempty(conc_unit{2}) % then the unit is part of the first string
    conc_unit = conc_unit{1}(21:end); % removes the text "Concentration unit: " to leave the unit itself
else % there is a tab between the text and the unit    
    conc_unit = conc_unit{2}; % second element is the unit
end

nxtline = fgetl(fileID); % reads next line from the file

if isempty(strfind(nxtline,'Concentration time')) && isempty(strfind(nxtline,'concentration time')) 
% This is not the recommended syntax for Matlab. However, <contains> was
% added in release 2016b, so this formulation works error free in older
% versions. Preferred syntax for newer versions:
% if ~contains(nxtline,'Concentration','IgnoreCase',true)
    error_str = sprintf('Format error in data set %1s: expected exposure scenario header below concentration unit',filename);
    return
end 
scen_header = strsplit(nxtline,'\t','CollapseDelimiters',0); % split text string at the tab and makes cell array
if length(scen_header) ~= length(surv_header)
    error_str = sprintf('Format error in data set %1s: exposure scenario must have same number of columns as survival data',filename);
    return
end

nxtline    = fgetl(fileID); % reads next line from file
scen_concs = []; % initialise an empty matrix for scenario information

while nxtline ~= -1 % continue until we find the end of file
    
    cell_line = strsplit(nxtline,'\t','CollapseDelimiters',0); % return cell array with expressions
    if length(scen_header) ~= length(cell_line)
        error_str = sprintf('Format error in data set %1s: exposure scenario does not have same number of columns at each time point',filename);
        return
    end
    
    scen_concs = [scen_concs; nan(1,length(cell_line))]; % initialise a line of NaNs
    for i = 1:length(cell_line) % run through all cells of the array
        scen_i = str2double(cell_line{i}); % convert string to number
        scen_concs(end,i) = scen_i; % note that any text fields will automatically become NaN
    end
    nxtline = fgetl(fileID); % reads next line from file
end 

fclose(fileID); % close the file

surv_header = surv_header(2:end); % remove the text in the first column
scen_header = scen_header(2:end); % remove the text in the first column

[tst,loc_scen] = ismember(surv_header,scen_header); % are all identifiers of survival also in the scenarios?
if any(tst==0)
    error_str = sprintf('Format error in data set %1s: some scenarios for survival do not have matching scenario identifier',filename);
    return
end

% If that check goes through, all survival identifiers have a matching
% scenario identifier. However, they may still be in the wrong order. Here,
% I decided to change the order of the scenario data set to match the
% survival data. However, it may be better to produce an error.

% scen_header = scen_header(loc_scen); % put scenario identifiers in same order as survival data
scen_concs(:,2:end) = scen_concs(:,loc_scen+1); % put scenario info in same order (first column is times)

% Collect all information into a structure <data> for output.
data.name     = surv_header;
data.surv     = surv_nmbrs;
data.scen     = scen_concs;
data.concunit = conc_unit;
