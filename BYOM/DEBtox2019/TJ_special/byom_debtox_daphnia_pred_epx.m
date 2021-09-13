%% BYOM, simple DEBtox model: byom_debtox_daphnia_pred_epx.m
%
% *Table of contents*

%% About
% * Author: Tjalling Jager
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* simple DEBtox model for toxicants, based on DEBkiss and
% formulated in compound parameters. The model includes flexible modules
% for toxicokinetics/damage dynamics and toxic effects. The DEBkiss e-book
% (see <http://www.debtox.info/book_debkiss.html>) provides a partial
% description of the model; a publication is in preparation that contains
% the full details. 
%
% *This script:* This script demonstrates the use of the parameter-space
% explorer from the openGUTS project (see <http://www.openguts.info/>) in
% making predictions for FOCUS profiles.
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn off the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(1) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% NOTE: time MUST be entered in DAYS for the estimation of starting to
% provide proper search ranges!

% scaled damage
DATA{1} = [0]; % there are never data for this state

% body length
DATA{2} = [0];

% cumulative reproduction (nr. offspring per mother)
DATA{3} = [0];

% survivors on each observation time
DATA{4} = [0];

%% Call the Matlab GUI open-file element to load profile and MAT file

[conf_type,fname_prof,par] = select_pred([2 1]); 

% Note that this function sets glo.mat_nm (to allow reading a sample from a
% different filename than the one indicated by the name of THIS script). It
% sets glo.Tbp as well. It sets it to the value used for the model fitting
% that generated the MAT file, or zero otherwise. When using an older MAT
% file (pre BYOM v6), make sure to set glo.Tbp to the correct value if it
% needs to be >0.

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [1]; % the scenarios, taken from the exposure profile (here: 1 as used for loaded files) 
X0mat(2,:) = 0; % initial values state 1 (scaled damage)
X0mat(3,:) = 0; % initial values state 2 (body length, initial value overwritten by L0)
X0mat(4,:) = 0; % initial values state 3 (cumulative reproduction)
X0mat(5,:) = 1; % initial values state 4 (survival probability)

% Put the position of the various states in globals, to make sure that the
% correct one is selected for extra things (e.g., for plotting in
% plot_tktd, in call_deri for accommodating 'no shrinking', for population
% growth rate).
glo.locD = 1; % location of scaled damage in the state variable list
glo.locL = 2; % location of body size in the state variable list
glo.locR = 3; % location of cumulative reproduction in the state variable list
glo.locS = 4; % location of survival probability in the state variable list

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
  
% global parameters as part of the structure glo
glo.FBV = 0.02;    % dry weight egg as fraction of structural body weight (-) (for losses with repro; approx. for Daphnia magna)
glo.KRV = 1;       % part. coeff. repro buffer and structure (kg/kg)
glo.kap = 0.8;     % approximation for kappa
glo.yP  = 0.8*0.8; % product of yVA and yAV (assume they are both 0.8)
glo.Lm_ref = 5;    % reference max length for scaling rate constants
glo.len = 1;       % switch to fit length 1) with shrinking, 2) without shrinking (used in call_deri.m)
% NOTE: the settings above are species specific! Make sure to use the same
% settings for validation and prediction as for calibration! 
% 
% NOTE: For arthropods, one would generally want to fit the model without
% shrinking (since the animals won't shrink in length). However, here it is
% turned off to allow you to identify if/when shrinking is triggered. 

% All parameters have been loaded from file using <select_pred>.

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% glo.t = linspace(0,Tend,100); % need to define time vector as we have no data

% specify the y-axis labels for each state variable
glo.ylab{1} = ['scaled damage (',char(181),'g/L)'];
glo.ylab{2} = 'body length (mm)';
if glo.Tbp > 0
    glo.ylab{3} = ['cumul. repro. (shift ',num2str(glo.Tbp),'d)'];
else
    glo.ylab{3} = 'cumul. repro. (no shift)';
end
glo.ylab{4} = 'survival fraction (-)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number
% Note: these legend labels will not be used when we make a glo.LabelTable

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 
% 
% NOTE: for this package, the options useode and eventson in glo will not
% be functional: the ODE solver is always used, and the events function as
% well.

glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances. Use 1 for quick analyses, but check with 3 to see if there is
% a difference! Especially for time-varying exposure, there can be large
% differences between the settings!
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% NOTE: for FOCUS profiles, breaking the time vector does not help accuracy
% much, but it does slow down the calculations a lot.

% No need for optimisation or standard plotting here.

%% Calculate ECx for all endpoints
% The function calc_ecx will calculate ECx values, based on the parameter
% values (and a sample from parameter space to make CIs). The ECx is
% calculated assuming constant exposure (which is in its definition) at the
% time points requested.

opt_conf.type    = conf_type; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 1; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_ecx.statsup  = []; % states to suppress from the calculations (e.g., glo.locS)
opt_ecx.par_read = 0; % when set to 1 read parameters from saved set, but do NOT make CIs

% calc_ecx(par,[1 5 10 15 21],opt_ecx,opt_conf);
% % Note that first input can be par_out, but if left empty it is read from
% % the saved file. Second input is a vector with time points on which ECx is
% % calculated (may extend beyond the time vector of the data set).

%% Analyse an exposure profile to check safety margins
% The functions below analyse the impacts of an exposure profile, given a
% certain exposure profile (e.g., from FOCUS). First, a 'moving
% time-window' analysis is done, with a fixed window of 21 days, and
% several fixed MFs. This will be a relatively fast and practical way to
% see if a certain critical safety margin is available. Next, the LPx can
% be calculated exactly.
% 
% Note that the functions calc_epx and calc_epx_window use the options
% structure opt_ecx. This option structure, by default, makes sure that the
% parameter called _hb_ is set to zero in the analysis (both in the
% best-fit parameter set and in the sample).

% opt_conf.type     = conf_type; % make intervals from 1) slice sampler, 2) likelihood region, 3) parspace explorer
opt_conf.type     = 0; % make intervals from 1) slice sampler, 2) likelihood region, 3) parspace explorer
opt_conf.lim_set  = 1; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
% opt_conf.n_lim    = 100; % size of limited set (likelihood-region and parspace only)
% opt_ecx.statsup   = [glo.locL]; % states to suppress from the calculations (e.g., glo.locL)
opt_ecx.rob_win   = 0; % set to 1 to use robust EPx calculation for moving time windows, rather than with fzero
opt_ecx.rob_rng   = [0.5 300 100 2]; % range for calculation of robust EPx, nr points, and log (0), normal (1), or custom scale (2)
% Note: when using the custom range (4-th element set to 2), the first 3
% elements of the vector will be ignored (settings in opt_ecx.rob_cust will
% be used).
% 
% Note: robust EPx calculation calculates EPx at the given steps only and
% then interpolates. This is not very precise, but with carefully selected
% steps it will be the lowest. In rather extreme cases, there will be more
% than one EPx for a specific window and a specific endpoint (when there
% are effects on growth, and feedbacks affecting k_d; a warning will be
% given). The 'regular' method may end up in either EPx (or produces an
% error), while 'robust' has a far better chance to yield the lowest.
% Robust is only used for EPx calculation, not for effect windows with
% fixed MFs.

opt_ecx.prune_win = 1; % set to 1 to prune the windows to keep the interesting ones
% Note: this setting screens all windows to find the window with highest
% minimum concentration. Any window whose maximum is lower than this value
% can be ignored. This is not guaranteed to work when there
% are effects on growth, and feedbacks affecting k_d (a warning will be
% given).
Twin = 21; % length of the time window (one element)

% -------------------------------------------------------------------------
% First, calculate effect (at end of each time window) for a range of fixed
% multiplication factors (MF).
opt_ecx.mf_range  = [3 10 30 100 300]; % range for fixed MFs to make plots with for calc_effect_window
calc_effect_window(par,fname_prof,Twin,opt_ecx,opt_conf);
% Note that first input can be par_out, but if left empty it is read from
% the saved file. This is needed now, as par_out is not known, since
% calc_optim is not called above.
% -------------------------------------------------------------------------
% We can also calculate explicit EPx values for each window. This will be
% considerably slower (especially when making CIs).
opt_ecx.Feff = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx/EPx
[MinColl,~,ind_traits] = calc_epx_window(par,fname_prof,Twin,opt_ecx,opt_conf);
% Note that first input can be par_out, but if left empty it is read from
% the saved file. This is needed now, as par_out is not known, since
% calc_optim is not called above.
% -------------------------------------------------------------------------

return

% Continue with a more detailed analysis of one window, including CI.
% MinColl contains all output, including most sensitive EPx and time point
% at which it occurs ... so we can automatically extract that time point.
[~,ind_min] = min(MinColl{opt_ecx.Feff==0.1}(:,2)); % find where lowest EP10 is 
Tstart = MinColl{opt_ecx.Feff==0.1}(ind_min,3); % start time for window with lowest EP10
% switch ind_traits(ind_min) % suppress states that are NOT the most sensitive ones
%     case glo.locS
%         opt_ecx.statsup = [glo.locL;glo.locR]; % states to suppress from the calculations (e.g., glo.locL)
%     case glo.locL
%         opt_ecx.statsup = [glo.locS;glo.locR]; % states to suppress from the calculations (e.g., glo.locL)
%     case glo.locR
%         opt_ecx.statsup = [glo.locL;glo.locS]; % states to suppress from the calculations (e.g., glo.locL)
% end

[~,ind_min] = min(MinColl{opt_ecx.Feff==0.5}(:,2)); % find where lowest EP10 is 
Tstart = MinColl{opt_ecx.Feff==0.5}(ind_min,3); % start time for window with lowest EP10

% Note: alternatively, you can set Tstart manually to the starting point of
% the window of choice 

opt_conf.type     = conf_type; % make intervals from 1) slice sampler, 2) likelihood region, 3) parspace explorer
opt_tktd.addzero  = 1; % set to 1 to always add a concentration zero to X0mat
opt_tktd.min      = 1; % set to 1 to show a dotted line for the control (lowest) treatment
opt_tktd.notitle  = 1; % set to 1 to suppress titles above plots
opt_tktd.max_exp  = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
% Note: for sub-lethal endpoints here, the dotted line for the control is
% handy. The overall title can be suppressed as it may block the titles of
% the sub-panels.

opt_ecx.batch_epx = 0; % when set to 1 use batch mode (no output to screen)
opt_ecx.Feff = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx/EPx

Trng = [Tstart Tstart+Twin]; % time range from the profile as specified in fname_prof to calculate EPx
% Note: the Trng has 2 elements, for the start and end time of the window

% Use one of the following functions to generate EPx. Here, the decision is
% based on the setting used for the moving-time window approach above.
if opt_ecx.rob_win ~= 1 % EPx with root finding
    calc_epx(par,fname_prof,Trng,opt_ecx,opt_conf,opt_tktd); 
else % EPx with stepwise increase in MF
    calc_epx_robust(par,fname_prof,Trng,opt_ecx,opt_conf,opt_tktd); 
end
% leave opt_tktd empty to NOT make a plot, or set opt_ecx.plot = 0
% leave Twin empty to use the entire profile as specified in fname_prof

