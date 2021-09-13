function disp_pmat(pmat,pmat_print,WRAP)

% Usage: disp_pmat(pmat,pmat_print,WRAP)
%
% This function prints <pmat> and <pmat_print> on screen in a formatted
% manner. This links to the parameter-space explorer algorithm, which
% prints the fitted parameters with their search ranges before start of the
% algorithm, and the best estimates with their CI after.
% 
% Inputs
% <pmat>        parameter matrix
% <pmat_print>  matrix with information on each parameter (best value and CI info)
% 
% Author     : Tjalling Jager 
% Date       : May 2020
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <disp_mat> code that is distributed
% as part of the Matlab version of openGUTS (see
% http://www.openguts.info). Therefore, this code is distributed under the
% same license as openGUTS (GPLv3). The modifications are only to ensure
% that the code operates in the general BYOM framework.
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

glo  = WRAP.glo;
glo2 = WRAP.glo2;

%% BLOCK 1. Initial things.
par_names = glo2.names; % names of the model parameters

n_par   = size(pmat,1);   % number of parameters in <pmat> in total
n_fit   = sum(pmat(:,2)); % number of fitted parameters
ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)

% put parameters that need to be fitted on log scale back on normal scale
pmat(pmat(:,5)==0,1) = 10.^(pmat(pmat(:,5)==0,1));
% NOTE THAT FIRST COLUMN OF PMAT IS HERE INPUT TO THIS FUNCTION, AND IS ON
% LOG-SCALE WHERE NEEDED (BOUNDS ARE STILL ON NORMAL SCALE).

%% Block 2. Print parameter values/ranges/intervals.

if isempty(pmat_print) % then we only print bounds (before an optimisation run)
    
    disp('=====================================================================')
    for i_p = 1:n_par % run through parameters
        if pmat(i_p,5) == 0 % create a string for log or normal scale fitting
            logsc = 'log';
        else
            logsc = 'norm';
        end
        if pmat(i_p,2) == 1 % for fitted parameters
            fprintf('%-5s bounds: %10.4g - %10.4g fit: %1.0f (%s)\n',par_names{i_p},pmat(i_p,3),pmat(i_p,4),pmat(i_p,2),logsc)
        else % for non-fitted parameters
            fprintf('%-5s fixed : %10.4g              fit: %1.0f\n',par_names{i_p},pmat(i_p,1),pmat(i_p,2))
        end
    end
    disp('=====================================================================')
    
else % then we print the results of a calibration run and more info needs to be printed
    
    % Replace flagged values in <pmat_print> by an ASCII code for display.
    prt_tst = 32*ones(n_fit,2); % make matrix with spaces (two columns)
    prt_tst(pmat_print(:,6:7)==1) = 42; % replace 'hitting bounds' flag by asterisk

    disp('==========================================================================')
    for i_p = 1:n_fit % run through FITTED parameters
        parnr = ind_fit(i_p); % index of this fitted parameter in <pmat>
        if pmat(parnr,5) == 0 % create a string for log or normal scale fitting
            logsc = 'log';
        else
            logsc = 'norm';
        end
        % Print the line to screen (broken up in three command lines).
        fprintf('%-5s best: %#10.4g (%#10.4g%c -',par_names{parnr},pmat_print(i_p,1),pmat_print(i_p,2),char(prt_tst(i_p,1)))
        fprintf('%#10.4g%c) ',pmat_print(i_p,3),char(prt_tst(i_p,2)))
        fprintf('fit: %1.0f (%s)\n',pmat(parnr,2),logsc)
    end
    disp('==========================================================================')
    if any(prt_tst(:) == 42) % we are printing at least one asterisk
        disp('* edge of 95% parameter CI has run into a boundary')
        disp('  (this may also have affected CIs of other parameters)')
    end
    if any(pmat_print(:,8) == 1) % there is at least one broken CI
        disp('==========================================================================')
        for i_p = 1:n_fit % run through parameters
            parnr = ind_fit(i_p); % index of this fitted parameter in <pmat>
            if pmat_print(i_p,8) == 1 % is there a flag for a broken CI?
                disp(['Confidence interval for parameter ',par_names{parnr},' is a broken set.'])
            end
        end
        disp('==========================================================================')
    end
    
end

%% Block 3. After calibration, we need to print some additional information

if ~isempty(pmat_print) % then we print the results of a calibration run and more info needs to be printed
    
    if isfield(glo,'loc_kd') % then we are probably doing a GUTS or DEBtox analysis
        
        loc_kd = glo.loc_kd; % location for kd
        loc_kd_fit = find(ind_fit==loc_kd); % location for kd in fitted pars
        
        if ~isempty(loc_kd_fit)
            disp(' ')
            disp('==========================================================================')
            % Calculate and print the DRT95 from <kd>.
            DRT95   = log(20)/pmat_print(loc_kd_fit,1); % calculate DRT95 from <kd>
            DRT95lo = log(20)/pmat_print(loc_kd_fit,3); % and its CI
            DRT95hi = log(20)/pmat_print(loc_kd_fit,2);
            fprintf('Depuration/repair time (DRT95)    : %#1.3g (%#1.3g - %#1.3g) days \n',DRT95,DRT95lo,DRT95hi)
            disp('Note that this depuration time does NOT account for any feedbacks from')
            disp('life history traits back to TK/damage dynamics (in DEBtox analyses).')
            disp('==========================================================================')
        end
    end
end

