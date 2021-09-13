%% BYOM, byom_logistic_sim.m, simulator example 
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
% *This script:* If you start building a model, simulations are useful to
% get a better idea of what the model can do. Make sure that the ODE
% version of the model is used (set _useode=1_).
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.
clear, clear global % clear the workspace and globals
global X0mat        % make scenarios global
global glo          % allow for global parameters in structure glo
diary off           % turn off the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine % set path to the BYOM/engine directory
glo.basenm = mfilename; % remember the filename for THIS file for the plots

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

% these scenarios specify different starting conditions: N(0)
X0mat(1,:) = [1 2 3 4 5];        % scenario identifiers
X0mat(2,:) = [10 12 50 100 150]; % initial value for state 1
   
%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = value
par.r   = 0.1;  % growth rate
par.K   = 100;  % carrying capacity
par.p   = 0;    % predation rate

%% What to calculate and what to plot
% Specify what to calculate and what to plot.

glo.t = linspace(0,100,500); % time vector for the model curves in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'population size (N)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (t)';
glo.leglab1 = 'Scen. '; % legend label before the 'scenario' number
glo.leglab2 = ''; % legend label after the 'scenario' number

% set options for call_deri here (we don't do prelim_checks for simulations)
glo.useode   = 1; % calculate model using ODE solver (1)
glo.eventson = 0; % events function on (1) or off (0)
glo.stiff    = 0; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)

%% Calculations and plotting
% Also try plottype 5 for this example! See the SIMbyom package for more
% detailed information on the simulator options.
 
opt_sim.plottype = 3; % select a plotting option
% 1) 3d plot, 
% 2) 2d plot, 
% 3) states in subplots, 
% 4) scenarios in subplots,
% 5) dx/dt versus x
% 6) 2d plot for each scenario separate

sim_and_plot(par,opt_sim); % call the script which calculates and plots (simulation only)

%% Local sensitivity analysis
% Local sensitivity analysis of the model. All model parameters are
% increased one-by-one by a small percentage. The sensitivity score is
% by default scaled (dX/X p/dp) or alternatively absolute (dX p/dp).
%
% Options for the sensitivity can be set using opt_sens (see
% prelim_checks.m).
 
prelim_checks % call prelim_checks to create some settings needed for the sensitivity analysis

% % UNCOMMENT LINE(S) TO CALCULATE
% opt_sens.type  = 1; % type of analysis 1) relative sensitivity 2) absolute sensitivity
% opt_sens.state = 0; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
% calc_localsens(par,opt_sens)