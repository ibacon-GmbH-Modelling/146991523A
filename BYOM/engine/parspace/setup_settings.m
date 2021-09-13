function SETTINGS_OPTIM = setup_settings(varargin)

% Usage: SETTINGS_OPTIM = setup_settings(varargin)
%
% This function is run several times when using the parameter-space
% explorer (parspace). It defines a number of settings for the
% parameter-space explorer. It would be better to have all these setting in
% the option structure <opt_optim>. However, there are quite a lot of them
% ... so for now, I will keep them in a separate function.
%
% Author     : Tjalling Jager
% Date       : February 2021
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <initial_setup> code that is
% distributed as part of the Matlab version of openGUTS (see
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

if numel(varargin)>0
    rough = varargin{1}; % rought option setting available as input argument
else
    rough = 1; % set to 1 to use rough settings, and 0 to use settings as in openGUTS
end

%% BLOCK 3. Define criteria for judging the likelihood ratio.
% Define the chi-square criteria for 1-5 df at 95% confidence (so we don't
% need <chi2inv> from the statistics toolbox). In Matlab this table is
% generated with <(chi2inv(0.95,1:6))'>.

% Critical values for 1-7 d.f. at alpha=0.05 (0.95 CI). Note: this table
% goes to df=7. However, for more than 5 parameters, the algorithm is
% unlikely to work properly!
SETTINGS_OPTIM.crit_table = [3.8415   
                  5.9915   
                  7.8147   
                  9.4877   
                 11.07    
                 12.592
                 14.0671];

SETTINGS_OPTIM.crit_table = min(SETTINGS_OPTIM.crit_table(5),SETTINGS_OPTIM.crit_table);
% TEST TEST this limits the volume of space that we are interested in (the
% open blue points in the parspace plot) to the critical value for df=5
             
% Next criterion selects the chi-square values to select the sub-sample to
% propagate to model predictions. In principle, this could be only the
% samples that sit on the criterium for 0.95, df=1. However, as we have a
% discrete (and limited) sample, it is best to take a bit more. I suggest
% taking a band around the value for 0.95 df=1. At this moment, it is set
% to the 92-96% area.

% SETTINGS_OPTIM.crit_prop = [3.0649 4.2179]; % Note: <[chi2inv(0.92,1) chi2inv(0.96,1)]>
        
SETTINGS_OPTIM.crit_prop = [-2*0.3764 0.3764] + SETTINGS_OPTIM.crit_table(1);
% NOTE: this has to be in line with the value for opt_conf.crit_add set in
% prelim_checks! (this is not done automatically, so take care) In
% <load_rnd>, 2 times the value is subtracted. This is very close to the
% openGUTS settings (it is now the 92.1-96% area).

%% BLOCK 5. Define settings for optimisation
% These settings are used in <calc_parspace>, and the functions that it
% calls. Note that not all tweakable options in <calc_parspace> are
% turned into a global setting (especially for the main round of sampling).
% I don't think there is much need to change those settings, and it would
% be hard to interpret their value outside of the context of the code.
% However, we could consider collecting them all as option settings in the
% future.

% Options for parameter optimisation (used in <calc_parspace>). Note that
% many of these settings are vectors: they change with each subsequent
% round of the analysis to allow the sample to contract.
SETTINGS_OPTIM.crit_add = [ 15  7  3   2   1.5   1  0.5  0.25]; % extra on top of chi2-criterion to select ok values for next round
SETTINGS_OPTIM.n_tr     = [NaN 60 40   30   20   10   10   10]; % number of random parameter tries in round 2 and further, for each ok parameter set
SETTINGS_OPTIM.f_d      = [NaN  1 0.65 0.42 0.27 0.18 0.12 0.08]; % maximum step for random mutations, as factor of initial grid spacing

% number of values that will at least continue to next round 
% NOTE: this is different from openGUTS as there always 200 is used. I here
% take this setting dependent on the number of fitted parameters.
% Otherwise, for 2 parameters (or 1) the initial round is pointless (all
% sets will continue to the next round). These values should be linked to
% the initial number of grid points!
SETTINGS_OPTIM.n_ok     = [5; % for 1 fitted parameters
   50 % for 2 fitted parameters
   200 % for 3 fitted parameters
   200 % for 4 fitted parameters
   400 % for 5 fitted parameters 
   800 % for 6 fitted parameters ...
   800]; % for 7 fitted parameters ...

% Stop criterion: minimum number of values within the total joint CI and
% the inner rim, which is made dependent on number of fitted parameters ...
SETTINGS_OPTIM.n_conf_all = [500  500   % for 1 fitted parameters
    2500 1500   % for 2 fitted parameters
    5000 2000   % for 3 fitted parameters
    7000 3000   % for 4 fitted parameters
    10000 5000  % for 5 fitted parameters
    15000 7500 % for 6 fitted parameters ...
    15000 7500]; % for 7 fitted parameters ...

SETTINGS_OPTIM.tries = 10; % initial grid points per fitted parameter
% Note: in openGUTS, the number of tries is between 8 and 14 depending on
% the parameter, but here we are in a general BYOM setting, so start with
% the same number for all parameters. In calc_parspace the setting in
% <tries> is copied over all fitted parameters. In theory, we could use a
% vector (different initial tries per parameter), but in BYOM, we don't
% know how many parameters there will be, and what they are ... We could
% use a setting in <glo> to optionally override <tries> with a vector ...

if rough >= 1% TEST TO MAKE IT FASTER (less points needed in cloud)
    SETTINGS_OPTIM.n_conf_all = [100  100   % for 1 fitted parameters
        250 150   % for 2 fitted parameters
        500 200   % for 3 fitted parameters
        700 300   % for 4 fitted parameters
        1000 500  % for 5 fitted parameters
        1500 750 % for 6 fitted parameters ...
        1500 750]; % for 7 fitted parameters ...
    SETTINGS_OPTIM.tries = 8; % initial grid points per fitted parameter
end
                        
% Settings for deciding on a restart for slow kinetics (used in <calc_parspace>).
SETTINGS_OPTIM.slowkin_corr = 0.70; % minimum correlation coefficient between <kd> and <mw>
SETTINGS_OPTIM.slowkin_pars = 0.05; % closeness of <kd> and <mw> to their lower bound (as fraction of total range)
% Settings for defining new upper bounds after restart (used in <calc_optim>)
SETTINGS_OPTIM.slowkin_f_mw = 3;  % factor by which to multiply current-highest <mw> to get new upper bound
SETTINGS_OPTIM.slowkin_f_kd = 10; % factor by which to multiply current-highest <kd> to get new upper bound

% Additional settings for the extra sampling (when profiling shows that
% sampling was poor in parts of the profile).
SETTINGS_OPTIM.crit_add_extra = 3; % continue with sets that are within this value from the MLL (this focusses on parameters within, or close to, the inner rim)
SETTINGS_OPTIM.f_d_extra      = 1; % initial maximum jump size as fraction of the grid spacing
SETTINGS_OPTIM.d_extra        = 0.1; % grid spacing as fraction of the parameter's total range (from joint 95% CI)
SETTINGS_OPTIM.gap_extra      = 0.25; % the gap distance between profile and sample that triggers resampling
SETTINGS_OPTIM.real_better    = 0.05; % difference between new optimum and old optimum (MLL) that leads to printing on screen

if rough >= 1% TEST TO MAKE IT FASTER (allow larger gaps)
    SETTINGS_OPTIM.gap_extra  = 2*0.25; % the gap distance between profile and sample that triggers resampling
end

% Settings for total number of mutation rounds (used in <calc_parspace>).
SETTINGS_OPTIM.n_max = 12;  % maximum number of rounds for the algorithm (default 12)
if rough > 1% TEST TO MAKE IT FASTER (allow larger gaps)
    SETTINGS_OPTIM.n_max = 4;  % maximum number of rounds for the algorithm (default 12)
    % This should be enough to find the best fit in most cases ...
end
