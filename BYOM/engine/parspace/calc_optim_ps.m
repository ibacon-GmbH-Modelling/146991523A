function [pmat,pmat_print,mll,names_out] = calc_optim_ps(pmat,opt_optim)

% Usage: [pmat,pmat_print,mll,names_out] = calc_optim_ps(pmat,opt_optim)
%
% Optimisation and finding a likelihood-based joint confidence region
% without the need for starting values. However, the parameter space to be
% searched must limited for the algorithm to be effective. Parameter space
% to be searched is defined as the min-max bounds of the parameters in the
% matrix <pmat>.
% 
% This function is a shell around the calibration algorithm in
% <calc_parspace>. It prepares for the optimisation, takes care of a
% restart (when slow kinetics is indicated during the optimisation), and it
% plots parameter space for saved sets if needed.
% 
% Note that the function <disp_pmat> is called to do the formatted printing
% (on screen) of <pmat> and <pmat_print>. The <pmat> is displayed before
% optimisation (to show the search ranges), the <pmat_print> after (with
% the optimised values and CIs).
%
% Inputs
% <pmat>        parameter matrix
% <opt_optim>   structure with options for optimisation
% 
% Outputs
% <pmat>        same as input matrix, but with best parameters in first column
% <pmat_print>  matrix with parameter info (best fit, CIs, etc.)
% <mll>         minus log-likelihood for the fit
% <names_out>   when using a saved set, this will return the saved names
%
% Author     : Tjalling Jager
% Date       : June 2021
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <calc_optim> code that is
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

global glo glo2

diary (glo.diary) % collect output in the diary "results.out"

SETTINGS_OPTIM = setup_settings(opt_optim.ps_rough); % read settings from file
% It would be better to have all these setting in the option structure
% <opt_optim>. However, there are quite a lot of them ... so for now, I
% will keep them in a separate function. The opt_optim.ps_rough allows for
% rough settings to be returned.

names_out = [];

%% BLOCK 1. Use the saved set, and don't make a new calculation.

if opt_optim.ps_saved == 1 % this is a flag for using a saved set
    
    % By default, this file will use the name of the base script to load a MAT
    % file. However, we also want to use saved MAT files for predictions in new
    % script files. The glo.mat_nm can be used for that.
    if isfield(glo,'mat_nm')
        filenm = glo.mat_nm;
    else
        filenm = glo.basenm;
    end
    
    pmat_in  = pmat; % remember input pmat
    glo2_rem = glo2;
    
    % Load the saved MAT file for this analysis (saved by <calc_parspace>).
    if exist([filenm,'_PS.mat'],'file') == 2 % check if it exists first.
        load([filenm,'_PS.mat'],'pmat','coll_all','pmat_print','coll_prof_pruned')
        % Load all of the saved information from the last <calc_parspace> run.
        
        % If we have parameter names in the saved set, we can use them for
        % plot_grid, rather than the names currently in glo2. This is
        % helpful for validation, as we might want to use a sample to make
        % predictions in a somehwat different model.
        variableInfo = who('-file',[filenm,'_PS.mat']);
        if ismember('names',variableInfo) % then we have a modern BYOM version that also saves the names of par
            load([filenm,'_PS.mat'],'names')
            if ~isequal(glo2.names,names)
                warning('Parameters from saved MAT file differ from those in the current workspace!')
                names_out = names;
            end
            glo2.names = names; % put names from SAVED pmat into glo2
        elseif size(pmat,1) ~= size(pmat_in,1)
            error('The parameter number in the saved set differs from that in the current workspace')
        end
        
    else % otherwise, produce an error
        error(['There is no confidence set with filename ',[filenm,'_PS.mat'],' found, so run calibration (with correct settings) first.'])
    end
    
    if isempty(coll_all) % it is empty if we saved a <pmat> with all parameters fixed (silly but possible)
        mll  = inf; % just return MLL as infinite
    else % otherwise, we have to reconstruct the the plot of parameter space from the saved sample
        mll  = coll_all(1,end); % extract maximum likelihood value from saved set (<coll_all> is sorted, so it is in the first row)
        plot_grid(pmat,[],coll_all,coll_prof_pruned,[],SETTINGS_OPTIM); % create plot but don't save plot again ... or should we?
    end
        
    % And display it to screen and diary.
    % n_fit = sum(pmat(:,2)); % number of fitted parameters
    disp('=================================================================================')
    disp('Displaying and plotting optimisation results from a saved set')
    disp('=================================================================================')
    disp('Best fit parameter set and 95% bounds from the saved set:')
    
    disp_pmat(pmat,pmat_print); % call dedicated function for printing the <pmat> and <pmat_print> (loaded from file)
    % NOTE THAT PMAT WAS SAVED ON LOG SCALE FOR LOG PARS!
    disp(' ')
    diary off   
    glo2 = glo2_rem;
    return % nothing further to do, so return to the main script

end

%% BLOCK 2. Call <calc_parspace> to do an optimisation.

% First, perform some checks on <pmat>
if any(pmat(:,3) == 0 & pmat(:,5) == 0)
    error('For log parameters, lower bound cannot be zero. Check the parameter structure.')
end
if any(pmat(:,4)-pmat(:,3)<=0 & pmat(:,2)==1)
    error('For at least one parameter, fitting is requested but the search range is of zero or negative width.')
end

if sum(pmat(:,2)) > 7
    error('For now, use this algorithm to fit 7 parameters or less!')
    % The function setup_settings is not prepared for more fitted
    % parameters, and the algorithm might not be stable.
elseif sum(pmat(:,2)) > 4
    warning('off','backtrace')
    warning('The parameter-space explorer is not thoroughly tested for more than 4 free parameters!')
    warning('on','backtrace')
end

disp(' ')
disp('Settings for parameter search ranges:')
disp_pmat(pmat,[]) % call <disp_pmat> for a formatted printing of the parameter search ranges
% For DEBtox analyses, it is nice to also show the moa and feedb settings,
% especially when running batch jobs.
if isfield(glo,'moa')
    disp(['Mode of action: ',num2str(glo.moa)])
end
if isfield(glo,'feedb')
    disp(['Feedbacks     : ',num2str(glo.feedb)])
end

tic % turn on the stopwatch to time the optimisations
% Call the optimisation function <calc_parspace> to do the optimisation.
[minmax,mll,stats,pmat_print,pmat] = calc_parspace(pmat,opt_optim,SETTINGS_OPTIM);
% NOTE: PMAT_PRINT IS ON NORMAL SCALE BUT PMAT ON LOG-SCALE WHERE NEEDED

% The <calc_parspace> will break the analysis if it runs into slow kinetics
% and when <mw> is estimated on normal scale. In that case, it needs to
% restart with <mw> on log-scale, and additionally, the parameter ranges
% can probably be decreased a bit.
if minmax(1) ~= -1 % then we stopped before the end ... this should only happen once (so no 'while' needed)
    ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
    loc_kd = glo.loc_kd; % location for kd
    loc_mw = glo.loc_mi; % location for mw (can be a vector if there are more thresholds)
    % Note: location of parameters in total parameter vector is defined in
    % <startgrid_guts> and <startgrid_debtox>.
    loc_kd_fit = find(ind_fit==loc_kd); % location for kd in fitted pars
    % loc_mw_fit = find(ind_fit==loc_mw); % location for mw in fitted pars
    loc_mw_fit = find(ismember(ind_fit,loc_mw)==1); % location for mw in fitted pars ADAPTED FOR POSSIBLY TWO THRESHOLDS (DEBTOX)
    SETTINGS_OPTIM = setup_settings(opt_optim.ps_rough); % read settings from file
        
    % NOTE: THIS NEEDS TO BE CHECKED FOR MORE THAN ONE THRESHOLD!
    loc_mw = loc_mw(pmat(loc_mw,2)==1); % only keep the ones that we fit in loc_mw!
    % this is needed as we'll modify pmat for these ones, and we might fit
    % only one of the thresholds. 
    
    disp(' ')
    % We suspect slow kinetics (currently the only reason for a restart).
    disp('Slow kinetics indicated: parameter-space explorer is restarting with these settings:')
    % <calc_parspace> returned the min-max of (relevant part of) parameter
    % cloud on normal scale.
    pmat(loc_mw,5) = 0; % for <mw>, put fit flag on log scale
    % Note: if someone *really* wants to keep <mw> on normal scale, there
    % is currently no simple way to avoid triggering a restart and a
    % log-scale for <mw>.
    
    % Note: for <kd>, we don't need to put it on log scale, as it already
    % is (see <calc_parspace>: if <kd> is on normal scale, there will be no
    % check on slow kinetics anyway, so no restart).
    
    % Adapt upper parameter bounds in <pmat> where possible, using the values in <minmax>.
    pmat(loc_mw,4) = min(pmat(loc_mw,4) , minmax(loc_mw_fit,2)*SETTINGS_OPTIM.slowkin_f_mw); % if possible, adapt upper boundary of <mw>
    pmat(loc_kd,4) = min(pmat(loc_kd,4) , minmax(loc_kd_fit,2)*SETTINGS_OPTIM.slowkin_f_kd); % also modify upper bound <kd> a bit if possible ...
    % NOTE: TRY TO MODIFY OTHER BOUNDS AS WELL?
    
    % Display new to-be-fitted parameter ranges on screen.
    disp_pmat(pmat,[]); % call dedicated function for printing the <pmat>
    
    % And restart the optimisation with <calc_parspace>.
    [~,mll,stats,pmat_print,pmat] = calc_parspace(pmat,opt_optim,SETTINGS_OPTIM);
end

%% BLOCK 3. Final printing on screen.
% BYOM will also print the standard output in <calc_optim>; this specific
% output (with CIs) will appear above that.

% The stats vector contains some interesting info from the analysis.
ind_final = stats(1); % index to the last element in the sample in the joint 95% CI
ind_inner = stats(2); % index to the last element in the sample in the inner rim
ind_prop  = stats(3); % number of elements in the sample that are used for propagating errors

disp(' ')
disp('=================================================================================')
disp('Results of the parameter-space exploration')
disp('=================================================================================')
disp(['   Sample: ',num2str(ind_final),' sets in joint CI and ',num2str(ind_inner),' in inner CI.'])
disp(['   Propagation set: ',num2str(ind_prop),' sets will be used for error propagation.'])
disp(['Best estimates and 95% CIs on fitted parameters'])
switch opt_optim.ps_rough
    case 1
        disp('Rough settings were used.')
        disp('  In almost all cases, this will be sufficient for reliable results.')
    case 2
        disp('VERY rough settings were used.')
        disp('  In most well-behaved cases, the best-fit parameter set will be reliable.')
        disp('  The sample will usually not be reliable for CIs.')
end
if opt_optim.ps_profs == 0 % then we skipped profiling and additional sampling!
    disp('Profiling and additional sampling steps were skipped.')
    disp('  CIs below are just the edges of the inner cloud (underestimating the true CIs).')
    disp('  Advised search ranges above are derived from approx. CIs +- 20%.')
end
disp_pmat(pmat,pmat_print); % call dedicated function for printing the <pmat>

disp(' ')
diary off  % close the diary that collected the output to screen
