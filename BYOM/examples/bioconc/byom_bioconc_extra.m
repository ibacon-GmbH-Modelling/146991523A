%% BYOM, byom_bioconc_extra.m, a detailed example 
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
% Fitting relies on the multinomial likelihood for survival and independent
% normal distributions (if needed after transformation) for continuous
% data. Results are shown on screen but also saved to a log file
% (results.out).
%
% In this example file (byom_bioconc_extra), we will walk step-by-step
% through a byom script file, explaining what happens and showing the
% outputs.
% 
% *The model:* An organism is exposed to a chemical in its surrounding
% medium. The animal accumulates the chemical according to standard
% one-compartment first-order kinetics, specified by a bioconcentration
% factor (_Piw_) and an elimination rate (_ke_). The chemical degrades at a
% certain rate (_kd_). When the external concentration reaches a certain
% concentration (_Ct_), degradation stops. This is useful to demonstrate the
% events function in call_deri.m, which catches this discontuity
% graciously. In the form of ODE's:
%
% $$ \frac{d}{dt}C_w=-k_d C_w \quad \textrm{as long as } C_w>C_t $$
%
% $$ \frac{d}{dt}C_i=k_e(P_{iw}C_w-C_i) $$
%
% *This script:* byom_bioconc_extra demonstrates fitting using replicated
% data, with different concentration vectors in both data sets.
% Furthermore, several options are explained, such as the inclusion of
% zero-variate data and priors for Bayesian analysis. More support in the
% BYOM manual downloadable from <http://www.debtox.info/byom.html here>. A
% stripped version (less text, less options) can be found as
% <byom_bioconc_start.html byom_bioconc_start.m>.
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Initial things
% Before we start, the memory is cleared, globals are defined, and the
% diary is turned off (output will be collected in the file results.out).
% The function pathdefine.m makes sure that the engine directory is added
% to the path. Make sure that this script is in a directory somewhere
% *below* the BYOM folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn off the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 to use parallel toolbox, if available)
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
%
% For each state variable there should be a data set in the correct
% position. So, _DATA{1}_ should be the data for state variable 1, etc. The
% curly braces indicate that this is a 'cell array'. If there are no data
% for a state, use _DATA{i}=0_ (if you forget, this will be automatically
% done in the engine script prelim_checks.m). Use NaN for missing data.
%
% Multiple data sets for each state can now be used. The first data set for
% state 1 becomes _DATA{1,1}_ and the second _DATA{2,1}_, etc. Note that
% each data set is treated as completely independent. For continuous data
% that means that each data set has its own error distribution (either
% treated as 'nuisance parameter' or provided by the user in _glo.var_).
% You can enter replicated data by adding columns with the same scenario
% number. The scenario number should occur only once in the initial values
% matrix _X0mat_.
%
% *Note:* for using survival data: enter the survival data as numbers of
% survivors, and not as survival probability. The model should be set up to
% calculate probabilities though. Therefore, the state variable is a
% probability, and hence the initial value in _X0mat_ below should be a
% probability too (generally 1). This deluxe script automatically
% translates the data into probabilities for plotting the results. For
% survival data, the weights matrix has a different meaning: it is used to
% specify the number of animals that went missing or were removed during
% the experiment (enter the number of missing/removed animals at the time
% they were last seen alive in the test).
%
% *Note:* In this example, the time vectors are not the same for both data 
% sets. This is no problem for the analysis. Also the concentration vectors
% differ, which is also accounted for, but take care to construct the
% matrix with initial states carefully to catch ALL scenarios.

% observed external concentrations (two series per concentration)
DATA{1} = [ 1  10 10 20 20 30 30
            0  11 12 20 19 28 29
            10  6  8 13 12 20 18
            20  4  5  9  7 10 12
            30  5  7  6  4  6  7
            40  4  5  6  5  4  3];
        
% weight factors (number of replicates per observation in DATA{1})
W{1} = [10 10 10 10 10 10
        10 10 10 10 10 10
         8  8  8  8  8  9
         8  7  7  8  8  9
         6  6  6  6  6  7];
    
% observed internal concentrations (two series per concentration)
DATA{2} = [0.5   10  10   20   20   50   50
           0    0   0    0    0    0    0
           10 570 600 1150 1100 2500 3000
           20 610 650 1300 1210 2700 2800
           25 590 620  960  900 2300 2200
           40 580 560  920  650 1200 1700];   

% if weight factors are not specified, ones are used for continuous data,
% and zero's for survival data.

%% Initial values for the state variables
% For each state variable, we need to specify initial values for each
% scenario that we want to fit or simulate in the matrix _X0mat_. The first
% row specifies the scenarios (here: exposure treatments) that we want to
% model. Here, there are 4 scenarios in _X0mat_, but each data set has only 3
% of them. Watch the plot to see how BYOM deals with that situation. 
% If you do not want to fit certain scenarios (exposure treatments) from
% the data, simply leave them out of _X0mat_.
%
% Note: if you do *not* want to start at _t_=0, specify the exact time in the
% global variable _glo.Tinit_ here (e.g., _glo.Tinit_ = 100;). If it is not
% specified, zero is used. The _X0mat_ thus defines the states at
% _glo.Tinit_, and not necessarily the value at the first data point.
% Plotting also always starts from _glo.Tinit_.
%
% Initial states, scenarios in columns, states in rows. First row are the
% identifiers of all scenarios (here: nominal concentrations). Second row
% is external concentration in mg/L, and third row body residues. For
% replicated data, scenarios should occur only once in _X0mat_.

X0mat = [10 20 30 50    % the scenarios (here nominal concentrations) 
          9 18 27 47    % initial values state 1 (actual external concentrations)
          0  0  0  0];  % initial values state 2 (internal concentrations)

% Note: The initial external concentrations are slightly different from 
% the "nominal" ones. 

% X0mat(:,3) = []; % for example, quickly remove the third scenario (conc 30)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. This means
% that you can address parameters by their name, instead of their position
% in a parameter vector. This makes the model definition in derivatives.m a
% lot easier. For each parameter, provide the initial value, whether you
% want to fit it or fix it to the initial value, the minimum bound, the
% maximum bound, and whether to fit the parameter on log10-scale or on
% normal scale. Fitting on log-scale is advisable for parameters that can
% span a very wide range; the optimisation routine can search this range
% more effectively if it is on log-scale. If you do not include this last
% value, a one will be filled in by prelim_checks.m (which means fitting on
% normal scale).

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.kd    = [0.1  1 1e-3 10  1];  % degradation rate constant, d-1
par.ke    = [0.2  1 1e-3 10  1];  % elimination rate constant, d-1
par.Piw   = [200  1 100  1e3 1];  % bioconcentration factor, L/kg
par.Ct    = [5    0 0    10  1];  % threshold external concentration where degradation stops, mg/L

% % Optionally, include an initial state as parameter (requires changes in call_deri.m too)
% par.Ci0   = [100  0 0 1e6];  % initial internal concentration (mg/L)

%% More options for the analysis
% Optionally, global parameters can be used in the structure glo (not
% fitted). However, see the text file reserved_globals.txt for names to
% avoid.

% % specify globals for parameters that are never fitted
% glo.Ct = 5; % the threshold could be made into a fixed global (also modify derivatives.m)
% % Special globals to modify the optimisation (used in transfer.m)
% glo.wts = [1 10];  % different weights for different data sets (same size as DATA)
% glo.var = [10 20]; % supply residual variance for each data set, after transformation (same size as DATA)

%% Zero-variate data and priors for Bayesian analyses
% Optionally, zero-variate data can be included. The corresponding model
% value needs to be calculated in call_deri.m, so modify that file too! You
% can make up your own parameter names in the structure _zvd_. Optionally,
% prior distributions can be specified for parameters, see the file
% calc_prior.m for the definition of the distributions. You must use the
% exact same names for the prior parameters as used in the _par_ structure.
% If you do not specify _zvd_ and/or _pri_ in your scripts, these options
% are simply not used in the analysis.
 
zvd.ku = [10 0.5]; % example zero-variate data point for uptake rate, with normal s.d.
 
% First element in pri is the choice of distribution.
pri.Piw   = [2 116 121 118]; % triangular with min, max and center
pri.ke    = [3 0.2 0.1];     % normal with mean and sd
% If no prior is defined, the min-max bounds will define a uniform one.
% Note that the prior is always defined on normal scale, also when the
% parameter will be fitted on log scale.

%% Time vector and labels for plots
% Specify what to plot. This also involves the time vector for the model
% lines (which may be taken differently than the time vector in the data
% set). If time vector glo.t is not specified, a default is constructed,
% based on the data set. You can specify the text that you would like to
% use for the axes, and the text used for the legend (which always includes
% the identifier for the scenario: the values in the first row of _X0mat_).

glo.t = linspace(0,50,100); % time vector for the model curves in days

% specify the y-axis labels for each state variable
glo.ylab{1} = 'external concentration (mg/L)';
glo.ylab{2} = 'internal concentration (mg/kg)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = 'mg/L'; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the functions are called that will do the calculation and the
% plotting. Note that calc_plot can provide all of the plotting information
% as output, so you can also make your own customised plots. This section,
% by default, makes a multiplot with all state variables (each in its own
% panel of the multiplot). When the model is fitted to data, output is
% provided on the screen (and added to the log-file results.out). The
% zero-variate data point is also plotted with its prediction (although the
% legend obscures it here).
%
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimisation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo.
%
% You can turn on the events function there too, to smoothly catch the
% discontinuity in the model. For the demo, the iterations were turned off
% (opt_optim.it = 0).

glo.eventson   = 1; % events function on (1) or off (0)
glo.useode     = 1; % calculate model using ODE solver (1) or analytical solution (0)
opt_optim.it   = 1; % show iterations of the optimisation (1, default) or not (0)
opt_plot.annot = 1; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend
opt_plot.bw    = 0; % if set to 1, plots in black and white with different plot symbols

% glo.diary = 'test_results.out'; % use a different name for the diary

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot);     % calculate model lines and plot them

% save_plot(gcf,'fit_example',[],3) % save active figure as PDF (see save_plot for more options)

%% Local sensitivity analysis
% Local sensitivity analysis of the model. All model parameters are
% increased one-by-one by a small percentage. The sensitivity score is
% by default scaled (dX/X p/dp) or alternatively absolute (dX p/dp).
%
% Options for the sensitivity can be set using opt_sense (see
% prelim_checks.m).
 
% % UNCOMMENT LINE(S) TO CALCULATE
% opt_sens.state = 2; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
% calc_localsens(par_out,opt_sens)

%% Profiling the likelihood
% By profiling you make robust confidence intervals for one or more of your
% parameters. Use the names of the parameters as they occurs in your
% parameter structure _par_ above. This can be a single string (e.g.,
% 'kd'), a cell array of strings (e.g., {'kd','ke'}), or 'all' to profile
% all fitted parameters. This example produces a profile for each parameter
% and provides the 95% confidence interval (on screen and indicated by the
% horizontal broken line in the plot).
%
% Notes: profiling for complex models is a slow process, so grab a coffee!
% If the profile finds a better solution, it breaks the analysis (as long
% as you keep the default opt_prof.brkprof=1) and displays the parameters
% for the new optimum on screen (and in results.out). For the NON-parallel
% version of calc_proflik, the new optimum is also immediately saves to the
% log file profiles_newopt.out. You can break of the analysis by pressing
% ctrl-c anytime, and use the values from the log file to restart (copy the
% better values into your _par_ structure). For parallel processing, saving
% would need more thought.
%
% Options for profiling can be set using opt_prof (see prelim_checks.m).

opt_prof.detail   = 1; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness
opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)

% % UNCOMMENT LINE(S) TO CALCULATE
% par_better = calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% if ~isempty(par_better)                 % if the profiling found a better optimum ...
%     calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
% end

%% Slice sampler
% The slice sampler can be used for a Bayesian analysis as it provides a
% sample from the posterior distribution. A .mat file is saved which
% contains the sample, to use later for e.g., intervals on model
% predictions. The output includes the Markov chain and marginal
% distributions for each fitted parameter. The function calc_conf can be
% used to put confidence intervals on the model lines. 
%
% Notes: MCMC is also slow process for complex models, so consider another
% coffee! If a prior is specified, it is plotted in the marginal posterior
% plot as a red distribution. Two types of credible intervals are
% generated:
%
% # Highest probability regions. This captures the 95% of the parameter
% values with the highest posterior probability. it is calculated by
% lowering the horizotal dotted line in the figures for each parameters,
% until the area between the cut-off points captures 95% of the density.
% This can lead to unequal probabilities in the tails.
% # Quantiles. The region between the 2.5 and 97.5% quantiles is taken.
% This leads to 2.5% probability in each tail. This is shown in the plots
% by the broken vertical lines.
% 
% Options for the slice sampling can be set using opt_slice (see
% prelim_checks.m). Options for confidence bounds on model curves can be
% set in opt_conf (see prelim_checks).
% 
% Taking on parameters on log scale helps, as the slice sampler seems quite
% bad at obtaining a representative sample when the parameters have very
% different ranges. Thinning will be needed to reduce the autocorrelation
% (which is too much here). Note the plot that is produced with details of
% the sample!

% % UNCOMMENT LINE(S) TO CALCULATE
% opt_slice.thin     = 10; % thinning of the sample (keep one in every 'thin' samples)
% opt_slice.burn     = 100; % number of burn-in samples (0 is no burn in)
% opt_slice.alllog   = 1; % set to 1 to put all parameters on log-scale before taking the sample
% calc_slice(par_out,100,opt_slice); % second argument number of samples (-1 to re-use saved sample from previous runs)

%% Likelihood region
% Another way to make intervals on model predictions is to use a sample of
% parameter sets taken from the joint likelihood-based conf. region. This
% is done by the function calc_likregion.m. It first does profiling of all
% fitted parameters to find the edges of the region. Then, Latin-Hypercube
% shooting, keeping only those parameter combinations that are not rejected
% at the 95% level in a lik.-rat. test. The inner rim will be used for CIs
% on forward predictions.
%
% Options for the likelihood region can be set using opt_likreg (see
% prelim_checks.m). For the profiling part, use the options in opt_prof.

opt_prof.detail  = 1; % detailed (1) or a coarse (2) calculation
opt_prof.subopt  = 0; % number of sub-optimisations to perform to increase robustness
opt_prof.re_fit  = 1; % set to 1 to automatically refit when a new optimum is found
opt_likreg.skipprof = 0; % skip profiling (and use profile from a saved mat file)

par_better = calc_likregion(par_out,500,opt_likreg,opt_prof,opt_optim); 
% Second entry is the number of accepted parameter sets to aim for. Use -1
% here to use a saved set.

if isstruct(par_better) % if the profiling found a better optimum ...
    print_par(par_better) % display on screen in formatted manner, so it can be copied into the code
    return % stop here, and don't go into plotting with CIs
end

%% Plot results with confidence intervals
% The following code can be used to make a standard plot (the same as for
% the fits), but with confidence intervals. Options for confidence bounds
% on model curves can be set using opt_conf (see prelim_checks).
% 
% Use opt_conf.type to tell calc_conf which sample to use: 
% -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% the model curves from the sample 
% 2) parameter sets from a joint likelihood region using the shooting 
% method (limited sets can be used), which will yield (asymptotically) 95% 
% CIs on predictions
% 3) as option 2, but using the parameter-space explorer

opt_conf.type    = 2; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time

out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs

% Here, we can also use the new plotting function for TKTD models. Even
% though this is not a TKTD model, we can still plot the internal
% concentration, with the treatments in separate panels.
glo.locC       = [1 2]; % tell plot_tktd that our first and second state variable are internal concentrations to plot
opt_tktd.repls = 0;     % plot individual replicates (1) or means (0)
opt_tktd.min   = 0;     % set to 1 to show a dotted line for the control (lowest) treatment

plot_tktd(par_out,opt_tktd,opt_conf);

%% Other files: derivatives
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file derivatives.m can be included in the
% published result.
% 
% <include>derivatives.m</include>

%% Other files: call_deri
% To archive analyses, publishing them with Matlab is convenient. To keep
% track of what was done, the file call_deri.m can be included in the
% published result.
%
% <include>call_deri.m</include>