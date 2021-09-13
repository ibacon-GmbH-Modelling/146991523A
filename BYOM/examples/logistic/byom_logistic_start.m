%% BYOM, byom_logistic_start.m, a quick fitting example 
%
% * Author: Tjalling Jager 
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* The standard logistic growth model, with parameters for the
% population growth rate (_r_) and the carrying capacity (_K_). Optionally,
% a constant predation rate (_p_) can be included. Model system:
%
% $$ \frac{dN}{dt} = r N \left(1-\frac{N}{K} \right) -p $$
%
% *This script:* Data set for algal growth from Schanz & Zahler (1981).
% Fitting with simplex optimisation and likelihood profiling to construct
% CIs. Note: you need to modify the initial parameter values and the
% intital states below to get a good fit!
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

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
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

% observed population size (in ml KMnO4)
DATA{1} = [1 1
    5	0.53
    7	1.29
    8	1.65
    9	2.25
    10	2.31
    11	2.91
    12	3.62
    13	4.11
    14	4.39
    16	4.67
    19	4.93];  

% if weight factors are not specified, ones are assumed in start_calc.m

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [1      % the scenario(s) to run
         1];  % initial values state 1: population size at t=0

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval];
par.r    = [0.1  1 0 5];  % population growth rate, d-1
par.K    = [10    1 0 20]; % carrying capacity, ml
% par.N0   = [1  1 0 1];  % initial population size, ml
par.p    = [0    0 0 1];  % predation rate, ml d-1

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'population size (ml KMnO4)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'Scen. '; % legend label before the 'scenario' number
glo.leglab2 = '';   % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the functions are called that will do the calculation and the
% plotting. 
% 
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. For the demo, the
% iterations were turned off (opt_optim.it = 0).

opt_optim.fit = 1; % fit the parameters (1), or don't (0)
opt_optim.it  = 1; % show iterations of the optimisation (1, default) or not (0)
glo.useode    = 1; % calculate model using ODE solver (1) or analytical solution (0)

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot);     % calculate model lines and plot them

%% Local sensitivity analysis
% Local sensitivity analysis of the model. All model parameters are
% increased one-by-one by a small percentage. The sensitivity score is
% by default scaled (dX/X p/dp) or alternatively absolute (dX p/dp).
%
% Options for the sensitivity can be set using opt_sense (see
% prelim_checks.m).
 
% % UNCOMMENT LINE(S) TO CALCULATE
% opt_sens.type  = 1; % type of analysis 1) relative sensitivity 2) absolute sensitivity
% opt_sens.state = 0; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
% calc_localsens(par_out,opt_sens)

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g., 'r'),
% a cell array of strings (e.g., {'r','K'}), or 'all' to profile all fitted
% parameters. This example produces a profile for each parameter and
% provides the 95% confidence interval (on screen and indicated by the
% horizontal broken line in the plot).
%
% *Note: for more post-calculations, see <byom_bioconc_extra.html
% byom_bioconc_extra.m>*
%
% Options for profiling can be set using opt_prof (see prelim_checks.m).

opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness

% % UNCOMMENT LINE(S) TO CALCULATE
% par_better = calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% 
% % Note: if a better optimum is found, the parameter structure is returned
% % in par_better. That structure could then be used in calc_and_plot, for
% % example.