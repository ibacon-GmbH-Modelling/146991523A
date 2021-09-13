function [conf_type,fname_prof,par,Tend,Cw] = select_pred(opt)

% Usage: [conf_type,fname_prof,par,Tend,Cw] = select_pred(opt)
% 
% This is a small helper function to allow the user to select MAT files for
% model curves with CIs, and exposure profiles for LPx/EPx predictions. 
% opt=1 to only load a MAT file
% opt=2 to also load an exposure profile. The exposure profile must be a 
%       txt file with time in the first column and concentration in the second.
% opt=[* 1] to also load the parameter structure par from file
% opt=[* 2] to also load the parameter structure par from file, and to use
%       the parameters for the second calibration data set (with glo.names_sep)
% opt = [* * 4] to also create an exposure profile with make_scen, using type 4
% 
% Outputs are:
% conf_type  the type of sample loaded (to be used to set opt_conf.type)
% fname_prof containing the filename for the exposure profile (to be used as input for the LPx/EPx predictions)
% par        parameter structure from saved set
% Tend       gives the last time point of the exposure profile (to define glo.t)
% Cw         the exposure profile itself (to be used with make_scen for test design)
% 
% Note that this function sets glo.mat_nm, which triggers load_rnd.m to
% load the selected mat file, rather than the one based on the name of the
% main script (in glo.basenm). Also note that this function sets glo.sel,
% glo.moa and glo.feedb, if these settings are included in the file name.
% This is used by some GUTS and DEBtox analyses. Note that this function
% sets glo.Tbp that specifes the brood-pouch delay used to generate the MAT
% file (only for specific DEBtox analyses). If par is requested as output,
% this function will also set glo.names_sep if it exists.
% 
% Author     : Tjalling Jager 
% Date       : August 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo

% use Matlab GUI to select MAT file for predictions
[fname_sample,filepath] = uigetfile('*.mat','Select a saved MAT file for predictions','MultiSelect','off'); 
if ~iscell(fname_sample) && numel(fname_sample) == 1 && fname_sample(1) == 0 % if cancel is pressed ...
    return % simply stop
end
if ~strcmp(filepath(1:end-1),pwd)
    error('The saved MAT file must be located in the present working directory')
end

if isfield(glo,'diary')
    diary (glo.diary) % collect output in the diary "results.out"
else
    diary ('results.out') % collect output in the diary "results.out"
end

disp('=================================================================================')
disp(['File selected for sample: ',fname_sample])
type = fname_sample(end-5:end-4); 
if strcmp(type,'MC')
    conf_type = 1; % Bayes slice sampler results
    disp('   Bayesian slice sample')
elseif strcmp(type,'LR')
    conf_type = 2; % likelihood region results
    disp('   likelihood-region sample')
elseif strcmp(type,'PS')
    conf_type = 3; % parspace explorer results
    disp('   parameter-space explorer sample')
elseif strcmp(type,'LP') 
    disp('   profile-likelihood results (does not include a sample!)')
    if exist([fname_sample(1:end-7),'_LR.mat'],'file') == 2
        conf_type = 2; % likelihood region results
        disp('   the similarly named likelihood-region file is available and can be used.')
    end
else
    error('The selected MAT file is not of known type.')
end

glo.mat_nm = fname_sample(1:end-7); % remove the last part (e.g., _PS.mat) of the filename so that the basename remains
% Note: setting this global implies that load_rnd will load this filename,
% rather than the one based on the script that generated it!

glo.Tbp = 0; % we must define glo.Tbp if it is not saved in the mat file!
% Note: this is only for DEBtox analyses ... but we'll keep it in for
% all analyses, as it does not hurt much.
listOfVariables = who('-file',fname_sample); % read stored variables in mat file
if ismember('Tbp', listOfVariables) % returns true if Tbp is stored
    load(fname_sample,'Tbp') % read glo.Tbp if it is saved in the mat file!
    glo.Tbp = Tbp; 
    disp(['   Loaded from saved MAT file, brood pouch delay  glo.Tbp = ',num2str(glo.Tbp)])
end

% Some GUTS analyses will save the death mechanism in the file name. The
% next section looks for the text '_sel' in the mat file name, and uses it
% to define glo.sel.
k = strfind(fname_sample,'_sel'); % find where _sel is in the filename
if ~isempty(k) % if it is not empty ...
    glo.sel = str2num(fname_sample(k+4)); % this is the death mechanism belonging to the saved file
    disp(['   Loaded from saved MAT file, death mechanism    glo.sel = ',num2str(glo.sel)])
end

% Some DEBtox analyses will save the MoA and the feedback configuration in
% the file name. This looks for the text '_moa' and '_feedb' to define
% glo.moa and glo.feedb.
ku = strfind(fname_sample,'_');
k1 = strfind(fname_sample,'_moa');
if ~isempty(k1) % if it is not empty ...
    moa_tmp = fname_sample(k1+4:ku(find(ku>k1,1,'first'))-1);
    % This reads from _moa up to the next underscore!
    for i = 1:length(moa_tmp)
        glo.moa(i) = str2num(moa_tmp(i));
    end
    disp(['   Loaded from saved MAT file, mode of action     glo.moa = ',num2str(glo.moa)])
end

k2 = strfind(fname_sample,'_feedb');
if ~isempty(k2) % if it is not empty ...
    feedb_tmp = fname_sample(k2+6:ku(find(ku>k2,1,'first'))-1);
    % This reads from _feedb up to the next underscore!
    for i = 1:length(feedb_tmp)
        glo.feedb(i) = str2num(feedb_tmp(i));
    end
    disp(['   Loaded from saved MAT file, feedback config. glo.feedb = ',num2str(glo.feedb)])
end    

if length(opt) > 1 && opt(2) > 0 % then we also need to return par!
    if ismember('par', listOfVariables) % returns true if par is stored
        load(fname_sample,'par') % read par if it is saved in the mat file!
        disp(['   Loaded from saved MAT file, parameter structure par'])
    else
        error('No parameter set is present in the saved MAT file!')
    end
    
    if ismember('names_sep', listOfVariables) % returns true if names_sep is stored
        load(fname_sample,'names_sep') % read glo.names_sep if it is saved in the mat file!
        glo.names_sep = names_sep;
        if ~isempty(names_sep)
            disp(['   Loaded from saved MAT file, extension namesep for ',num2str(length(glo.names_sep)),' parameters'])
        end
        if opt(2) > 1 % if the request is for a data set that is not the first ...
            disp(['   Parameters used in par are for calibration data set: ',num2str(opt(2))])
            % use the parameters for data set in opt(2)
            for i_sep = 1:length(names_sep) % run through extra parameter names for separate sets
                eval(['par.',names_sep{i_sep},' = par.',[names_sep{i_sep},num2str(opt(2)-1)],';']); % copy extra parameter level to par
            end
        end
    end
else % make sure the outputs are defined
    par = [];
end

if opt(1) > 1 % we also need to load an exposure profile
    
    % use Matlab GUI to load exposure profile from text file for predictions
    [fname_prof,filepath] = uigetfile('*.txt','Select a text file with exposure profile for predictions','MultiSelect','off');
    if ~iscell(fname_prof) && numel(fname_prof) == 1 && fname_prof(1) == 0 % if cancel is pressed ...
        return % simply stop
    end
    
    disp(['File selected for exposure profile: ',fname_prof])
    % if the profile is not located in the working directory, return the entire path
    if ~strcmp(filepath(1:end-1),pwd)
        disp(['   in directory: ',filepath])
        fname_prof = [filepath,fname_prof];
    end
    
    if isfield(glo,'int_scen') % scenarios have already been defined
        make_scen(-5,-1); % remove all spline info for exposure profiles (just to be on the safe side)
        disp('   previously defined exposure scenarios are deleted!')
    end
    
    % Load exposure profile from file to use with linear interpolation.
    % This is needed here to define a time vector glo.t.
    Cw   = load(fname_prof);
    disp(['   exposure matrix with ',num2str(size(Cw,1)),' time points and ',num2str(size(Cw,2)-1), ' scenarios'])
    Tend = Cw(end,1); % end of the profile (for defining glo.t)
    if size(Cw,2) == 2 && Cw(2,1) ~= 0 % if the second time element is a zero, we have a file defined with identifiers
        Cw = [1 1;Cw]; % this add a first row with scenario identifier
        disp('   the exposure scenario gets the identifier 1')
    end
    
    glo.t = linspace(0,Tend,max(100,2*size(Cw,1))); % also predefine glo.t
    disp(['   time vector glo.t predefined as 0-',num2str(Tend),' with ',num2str(length(glo.t)),' points'])
    if length(opt) > 2 % then we'll also define the exposure scenario here
        make_scen(opt(3),Cw); % create exposure profile with opt(3) as type
        disp('   exposure scenario defined in global (no need to call make_scen again)')
    end
    
else % make sure the outputs are defined
    fname_prof = [];
    Cw         = [];
    Tend       = [];
end


disp('=================================================================================')
diary off  % close the diary that collected the output to screen
