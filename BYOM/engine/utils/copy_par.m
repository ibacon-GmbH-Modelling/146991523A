function par = copy_par(par,par_saved,varargin)

% This copies the fitted parameters from a saved set into the par structure
% that will be used for plotting. It also makes sure that the fit marks are
% taken from the saved set. This is handy for validation in the debtox2019
% package: the control data from the validation set can be fitted, and used
% with the tox parameters from calibration (taken from a saved MAT file).
% There may be parameter differences between both sets since the debtox2019
% package allows for parameters to differ between data sets (which implies
% extra parameters).
%
% The varargin can be used to ask for specific non-fitted parameters to be
% copied from the saved file. This is useful when tox parameters where
% fixed in calibration, and need to be used for validation. That way, the
% user does not have to remember which fixed value was used!
% 
% Author     : Tjalling Jager 
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

names        = fieldnames(par); % extract all field names of par
names_saved  = fieldnames(par_saved); % extract all field names of par_saved from MAT file
ind_fittag   = ~strcmp(names_saved,'tag_fitted');
names_saved  = names_saved(ind_fittag); % make sure that the fit tag is not in names_saved
ind_fittag   = ~strcmp(names,'tag_fitted');
names        = names(ind_fittag); % make sure that the fit tag is not in names_saved

skip_fitmark = 0;
names_extra = {};
if ~isempty(varargin)
    if ischar(varargin{1}) % when entering extra strings ...
        names_extra = varargin; % varargin contains extra non-fitted parameters that need copying
    elseif isnumeric(varargin{1})
        skip_fitmark = varargin{1}; % when 1, don't put fitmark to zero when copying non-fitted parameters
    end
end

for i = 1:length(names_saved) % run through all parameters of *saved* par
    if ismember(names_saved{i},names) % the saved parameters is also in *input* par
        if eval(['par_saved.',names_saved{i},'(2) == 1;']) % the saved parameter is fitted
            eval(['par.',names_saved{i},' = par_saved.',names_saved{i},';']); % copy parameter from saved par to par
            % This copies the entire row, including fit mark, range, and log-setting
        else % the saved parameter is not fitted
            if skip_fitmark == 0
                eval(['par.',names_saved{i},'(2) = 0;']); % make sure fit mark in *input* par is also off
            end
            if ismember(names_saved{i},names_extra) % then it is in the additional parameters to be copied
                eval(['par.',names_saved{i},' = par_saved.',names_saved{i},';']); % copy parameters from saved par to par
            end
            
        end
    end
end