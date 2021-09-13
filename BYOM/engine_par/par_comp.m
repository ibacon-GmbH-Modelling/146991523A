function par_comp(par,par_tmp,varargin)

% This function compares the _par_tmp_ that it is entered into a certain
% function with the _par_ that is in a saved set.
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2

WRAP.glo  = glo;
WRAP.glo2 = glo2;

names      = fieldnames(par); % extract all field names of par (global)
ind_fittag = ~strcmp(names,'tag_fitted');
names      = names(ind_fittag); % make sure that the fit tag is not in names_tmp
nfields    = length(names);

names_tmp   = fieldnames(par_tmp); % extract all field names of par (global)
ind_fittag  = ~strcmp(names_tmp,'tag_fitted');
names_tmp   = names_tmp(ind_fittag); % make sure that the fit tag is not in names_tmp
nfields_tmp = length(names_tmp);

pmat_tmp = packunpack(1,par_tmp,0,WRAP); % transform structure *from input* into a regular matrix
pmat     = packunpack(1,par,0,WRAP); % transform structure *from saved set* into a regular matrix

warning('off','backtrace')
if nfields ~= nfields_tmp
    warning('The parameter matrix in the workspace is of a different size than in the saved set!')
    warning('off','backtrace'), disp(' ')
    return
end

% Identify parameter(s) to include in the comparison
ind_comp = true(size(pmat,1),1); % no parameters to ignore, unless something was entered
if ~isempty(varargin) 
    set_zero = varargin{1}; % then it is a cell array with the parameter(s) to set to zero
    if ~isempty(set_zero) % we may want to make a parameter zero (esp. background mortality)        
        if ~iscell(set_zero) % for backward compatibility
            set_zero = {set_zero}; % turn it into a cell array with one element
        end
        loc_zero = false(length(names),1); % create logical index with zeros
        for i = 1:length(set_zero)
            loc_zero = loc_zero == 1 | strcmp(names,set_zero{i})==1; % add this parameter to loc_zero (logical indexing)
        end
        ind_comp = ~loc_zero; % parameters to use for the comparison
    end
end

prpar = 0; % don't print the parameter set on screen, unless there is reason to
if ~isequal(pmat(ind_comp,:),pmat_tmp(ind_comp,:))
    
    disp(' ')
    if ~isequal(pmat(ind_comp,5),pmat_tmp(ind_comp,5)) % ah, it's just the log settings
        warning('Log settings of parameters in the saved set differ from that in the workspace.')
        prpar = 1;
    end
    if ~isequal(pmat(ind_comp,3:4),pmat_tmp(ind_comp,3:4)) % ah, it's just the boundaries
        warning('Boundaries of parameters in the saved set differ from that in the workspace.')
        prpar = 1;
    end
    if prpar == 1
        disp('         (differences in log-settings or boundaries can be expected in various analyses)')
    end
    if ~isequal(pmat(ind_comp,2),pmat_tmp(ind_comp,2))
        warning('Selection of fitted parameters is different in the saved set.')
        prpar = 1;
    end
    if ~isequal(pmat(ind_comp,1),pmat_tmp(ind_comp,1))
        % perhaps add some code here to NOT trigger this warning when the
        % difference is tiny (rounding errors); but take care of parameters
        % that may be zero ...
        
        diff_p  = abs(pmat(ind_comp,1)-pmat_tmp(ind_comp,1))./pmat(ind_comp,1); % relative difference
        diff_p2 = abs(pmat(ind_comp,1)-pmat_tmp(ind_comp,1)); % absolute difference
        ind_zero = find(pmat(ind_comp,1)==0 | pmat_tmp(ind_comp,1)==0); % any zero in either matrix?
        diff_p(ind_zero) = diff_p2(ind_zero); % use absolute difference instead
        
        if any(diff_p > 1e-5)
            warning('Best-fitting parameters differ from the ones in the workspace. The one in the workspace is used for the best curve, and the set from the saved set is used for the confidence intervals!')
            prpar = 1;
        else
            warning('Best-fitting parameters differ from the ones in the workspace, but the difference is negligible.')
        end
    end
    
    if prpar == 1
        fprintf('Parameter values from saved set, for checking \n');
        fprintf('=================================================================================\n');
        
        for i = 1:nfields % display all parameters on screen
            if pmat(i,5) == 0
                fprintf('%-6s %10.4g (fit: %1.0f) bounds: %6.4g - %6.4g log-scale \n',names{i},pmat(i,1),pmat(i,2),pmat(i,3),pmat(i,4))
            else
                fprintf('%-6s %10.4g (fit: %1.0f) bounds: %6.4g - %6.4g \n',names{i},pmat(i,1),pmat(i,2),pmat(i,3),pmat(i,4))
            end
        end
        fprintf('=================================================================================\n');
        fprintf('  \n');
    end
end
warning('off','backtrace'), disp(' ')