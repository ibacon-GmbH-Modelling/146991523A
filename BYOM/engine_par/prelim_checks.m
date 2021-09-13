% This script defines a number of useful variables, and does some
% preliminary checks on the data set. This script is called from the main
% script directly.
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo2

rng('shuffle')     % make sure that a *new* random seed is taken when restarting!
% rng(7)             % make sure that the *same* random seed is taken when restarting!
kill_waitbars = 1; % if set to 1, prelim_checks will kill all open waitbars (this needs to be turned off for some analyses I am running)

% This checks whether you have the parallel toolbox installed. If it is
% installed, it calculates how many cores to use for parallel computation.
% Note that <parfor> loops will be executed as <for> loops without this
% toolbox.
n_cores = 0; % by default, we have zero cores
% this is to avoid errors when running this code without the toolbox
v = ver; % take the version information
if any(strcmp('Parallel Computing Toolbox',{v.Name}))
    % If a parallel pool has not yet been created, create one with n_cores-1 workers
    n_cores = feature('numcores');
    if n_cores < 7 % number of physical cores available
        n_cores = n_cores-1; % this leaves 1 core for other work
    else % if you have at least 8 cores, leave 2 cores free
        n_cores = n_cores-2; % this leaves 2 cores for other work
        % on my 10-core workstation, using 9 does not improve speed over
        % using 8 cores. Using 8 beats using 6 cores, though.
    end
    % Note: the pool is not automatically deleted after the analysis, but
    % closes after 30 min idle, by default. Otherwise, kill it manually:
      % poolobj = gcp('nocreate'); % If no pool, do not create new one.
      % delete(poolobj) % force parallel pool to end
end
glo2.n_cores = n_cores; % put the number of cores to use in glo2

%% Defines options structures
% First, the options structure for the various calculation routines are
% specified with default values. 

% Options for parameter optimisation (used in calc_optim)
opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
opt_optim.it       = 1; % show iterations of the simplex optimisation (1, default) or not (0)
opt_optim.type     = 1; % optimisation method 1) simplex, 2) simulated annealing (experimental), 
                        % 3) swarm (experimental), 4) parameter-space explorer
opt_optim.simno    = 2; % for simplex: number of runs of fminsearch (starting again from previous best)
opt_optim.anntemp  = 1e-4; % for simulated annealing, stopping temperature
opt_optim.swno     = 200; % for swarm optimisation, number of particles
opt_optim.swit     = 20; % for swarm optimisation, number of iterations
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
opt_optim.ps_plots = 1; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
opt_optim.ps_notitle = 0; % set to 1 to suppress plotting title on parameter-space plot
opt_optim.ps_slice = 0; % set to 1 to use slice sampler for main rounds (for use by Tjalling only!)

% Options for plotting (used in calc_and_plot)
opt_plot.zvd     = 1; % turn on the plotting of zero-variate data (if defined)
opt_plot.sub     = 1; % switch for putting state variables as sub-plots into a figure (set to 1)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.cn      = 1; % if set to 1, connect model line to points (only for bw=1)
opt_plot.limax   = 0; % if set to 1, limit axes to the data set for each stage
opt_plot.sho     = 1; % set to 1 to show all scenarios, 0 to only plot model for scenarios with data
opt_plot.repls   = 1; % set to 1 to plot replicates, 0 to plot mean responses, 2 to make boxplots
opt_plot.notitle = 0; % set to 1 to suppress plotting titles on graphs with fits
opt_plot.outl    = 1; % set to 1 to identify outliers in the plot (points with weight zero; not for survival)
opt_plot.legsup  = 0; % set to 1 to suppress legends on fits
opt_plot.annot   = 0; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend
opt_plot.simulik = 1; % if set to 1, calculate the log-likelihood when fit=0 (in simulations)
opt_plot.statsup = []; % vector with states to suppress in plotting fits
opt_plot.transf  = 1;  % set to 1 to calculate means and SEs including transformations (if repls=0)
opt_plot.y_zero  = 1;  % set to 1 to force y-axis to start at zero

% Options for local sensitivity analyses (used in local_sens)
opt_sens.type     = 1; % type of analysis 1) relative sensitivity 2) absolute sensitivity
opt_sens.state    = 0; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
opt_sens.step     = 0.05; % fraction change in each parameter's value
opt_sens.notitle  = 0; % set to 1 to suppress plotting titles on graphs

% Options for asymptotical standard errors (used in calc_ase)
opt_ase.step      = 0.01; % relative stepsize

% Options for profiling likelihoods (used in calc_proflik)
opt_prof.detail   = 2; % detailed (1) or a coarse (2) calculation
opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness
opt_prof.subrng   = 5; % maximum factor on parameters (both higher and lower) for sub-optimisations
opt_prof.brkprof  = 0; % when a better optimum is located, stop (1) or automatically refit (2)
opt_prof.verbose  = 1; % set to 0 to suppress output to screen, 2 to only suppress plots
opt_prof.subann   = 0; % set to 1 to use simulated annealing (followed by simplex) instead of suboptimisations
                       % setting this options means that subopt has no effect
opt_prof.saved    = 0; % set to 1 to plot/display saved profiles

% Options for the slice sampler (used in calc_slice)
opt_slice.thin     = 1; % thinning of the sample (keep one in every 'thin' samples)
opt_slice.burn     = 100; % number of burn-in samples (0 is no burn in)
opt_slice.slwidth  = 10; % initial width of the slice (Matlab default is 10)
opt_slice.alllog   = 0; % set to 1 to put all parameters on log-scale before taking the sample
opt_slice.testing  = 1; % make additional tests on the sample: moving average and autocorrelation
opt_slice.subplt   = 1; % create single plot with subplots (1) or many individual plots (0)

% Options for the calculation of confidence intervals (used in calc_conf and plot_guts)
opt_conf.type     = 1; % use values from slice sampler (1), likelihood region (2), or parspace explorer (3) to make intervals
opt_conf.samerr   = 0; % include sampling error in bounds for survival data (set samerr=1; requires statistics toolbox, Bayes only)
opt_conf.n_samerr = 20; % number of sub-sampling trials for each parameter set (if samerr=1)
opt_conf.sens     = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
opt_conf.state    = 0; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
opt_conf.lim_set  = 0; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_conf.n_lim    = 200; % size of limited set (likelihood-region and parspace only)
opt_conf.crit_add = 0.3764; % small addition to chi2 criterion for inner rim (to make sure coverage is adequate)
                         % 0.3764 is in line with the value used in the openGUTS software
opt_conf.set_zero = {}; % parameter name (as cell array of strings) to set to zero for calc_conf (e.g., the background hazard)
opt_conf.use_par_out = 0; % set to 1 to use par as entered into the plotting function for CIs, rather than from saved set

% Options for the calculation of LCx/LPx values (used in calc_lcx_lim and calc_lpx_lim)
opt_lcx_lim.Feff  = 0.50; % effect level (>0 en <1), x/100 in LCx/LPx
opt_lcx_lim.plot  = 1; % set to 0 to NOT make a plot of LCx vs time or survival vs time at LPx
opt_lcx_lim.scen_plot = 1; % set to 0 to NOT make an extra plot of the exposure profile
opt_lcx_lim.scen_type = 4; % type of definition for the exposure scenario (as used in make_scen), if called with a file name
opt_lcx_lim.notitle   = 0; % set to 1 to suppress titles above ECx plots

% Options for the calculation of ECx and EPx values (used in calc_ecx, calc_epx and calc_epx_window)
opt_ecx.Feff      = [0.10 0.50]; % effect levels (>0 en <1), x/100 in ECx/EPx
opt_ecx.plot      = 1; % set to 0 to NOT make a plot of ECx vs time (or effects at different LPx)
opt_ecx.backhaz   = 'hb'; % parameter name (as string) to set to zero for LCx/LPx calculation to remove background mortality
opt_ecx.setzero   = {}; % parameter names (as cell array of strings) for extra parameters to be set to zero
opt_ecx.statsup   = []; % states to suppress from the calculations (e.g., locS)
opt_ecx.mf_range  = [1 10 100 1000]; % range for MFs to make plots with for calc_epx_window
opt_ecx.par_read  = 0; % when set to 1 read parameters from saved set, but do NOT make CIs
opt_ecx.batch_epx = 0; % when set to 1 use batch mode (no output to screen)
opt_ecx.notitle   = 0; % set to 1 to suppress titles above ECx plots
opt_ecx.start_neg = 1; % set to 1 to start the moving window at minus window width
opt_ecx.Cminmax   = [1e-7 1e6]; % absolute boundaries for ECx (and its CI), [min>0,max]
opt_ecx.calc_rgr  = 0; % set to 1 to calculate rgr in calc_effect_window (experimental!)
opt_ecx.rob_win   = 0; % set to 1 to use robust EPx calculation for moving time windows, rather than with fzero
opt_ecx.rob_rng   = [0.5 300 300 1]; % range for calculation of robust EPx, nr points, and log (0), normal (1), or custom scale (2)
opt_ecx.rob_cust  = [0.1:0.1:20 21:1:100 100:3:300]'; % custom range for ERA purposes (smaller steps in interesting region)
opt_ecx.nomarker  = 0; % set to 1 to suppress markers for the ECx-time plot
opt_ecx.mf_crit   = 30; % MF trigger for flagging potentially critical profiles (EP10<mf_crit)
opt_ecx.prune_win = 0; % set to 1 to prune the windows to keep the interesting ones

% Options for the joint likelihood confidence region (used in calc_likregion)
opt_likreg.skipprof = 0; % skip profiling step; use boundaries from saved likreg set (1) or profiling (2)
opt_likreg.chull    = 0; % set to 1 to plot convex hull that approximates 95% edges of the likelihood region
opt_likreg.axbnds   = 1; % bind axes on the bounds of the hyperbox (1), accepted sample (2), or inner region (3)
opt_likreg.burst    = 100; % number of random samples from parameter space taken every iteration
opt_likreg.lim_out  = 0; % set to 1 to sample from a smaller part of space (enough for forward predictions)
                         
% Options for plotting survival results for the GUTS packages (used in plot_guts) 
opt_guts.timeresp  = 1; % set 1 for multiplot with survival versus time for each treatment
opt_guts.doseresp  = 1; % set 1 for multiplot with survival versus concentration for each time point
opt_guts.timedeath = 1; % set 1 for multiplot with deaths-per-interval versus time for each treatment
opt_guts.single_dr = 0; % set 1 for single plot with survival versus concentration for each time point
opt_guts.single_td = 0; % set 1 for single plot with expected and observed deaths in each time interval
opt_guts.bw        = 0; % plot the single dose-response plot in black and white
opt_guts.cn        = 1; % use connecting lines between data and model curve when opt_guts.bw=1

% Options for the calculation of the intrinsic population growth rate
opt_pop.fscen   = [1 0.9 0.8]; % three scenarios with limited food
opt_pop.trep    = [0]; % vector with repro events (first one must be zero, only zero is continuous repro)
opt_pop.plt_fly = 1; % set to 1 for plotting on the fly

% Options for TKTD plotting
opt_tktd.repls    = 0; % plot individual replicates (1) or means (0)
opt_tktd.obspred  = 1; % plot predicted-observed plots (1) or not (0), (2) makes 1 plot for multiple data sets
opt_tktd.preds    = 0; % set to 1 to only plot predictions from X0mat without data
opt_tktd.addzero  = 0; % set to 1 to always add a concentration zero to X0mat
opt_tktd.max_exp  = 1; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.notitle  = 0; % set to 1 to suppress titles above plots
opt_tktd.transf   = 1; % set to 1 to calculate means and SEs including transformations
opt_tktd.min      = 1; % set to 1 to show a dotted line for the control (lowest) treatment
opt_tktd.statsup  = []; % states to suppress from the plots (e.g., locS)
opt_tktd.flip     = 0; % set to 1 to flip row 1 and 2 of the plot around
opt_tktd.plotexp  = 1; % set to 1 to plot exposure profile as area in damage plots
opt_tktd.lim_data = 0; % set to 1 to limit axes to data
opt_tktd.set_ctrl = 0; % set to 1 to use separate control per data set for the dotted lines
opt_tktd.set_no   = 1; % set to 1 to specify 'Set i' in the titles of the subplots
opt_tktd.sppe     = 0; % set to 1 to calculate SPPEs (relative error at end of test)
opt_tktd.symsz    = 6; % size for the plotting symbol (default 6)
opt_tktd.pltbar   = [];  % plot a bar at specified location, for selected scenarios

tic % turn on the stopwatch

if kill_waitbars == 1
    % Look if there are any waiting bars left open (might happen if a run did
    % not terminate properly).
    F = findall(0,'type','figure','tag','TMWWaitbar');
    delete(F) % close them
end

%% Definition helpful globals

if ~isfield(glo,'diary') % new option to modify the name of the diary per run
    glo.diary = 'results.out'; % by default, filename is "results.out"
end
if ~exist('pri','var') || isempty(pri)
    glo2.namesp = [];
    glo2.pri    = [];
else
    glo2.namesp = fieldnames(pri);  % extract all field names of prior (global)
    glo2.pri    = pri; % copy the prior definition into glo!
end
if ~exist('zvd','var') || isempty(zvd)
    glo2.namesz = [];
    glo.zvd     = [];
else
    glo2.namesz = fieldnames(zvd);  % extract all field names of zero-var data (global)
    glo.zvd     = zvd; % copy the zero-variate data set into glo!
end
if ~isfield(glo,'names_sep')
    glo.names_sep = {}; % start with empty global for names of data-set-specific parameters 
end

% In some cases, DATA and W may not have been defined (e.g., using
% simulations), but we may still need this script to allow for sensitivity
% analysis.
if ~exist('DATA','var')
    DATA = [];
end
if ~exist('W','var')
    W = [];
end

[n_D,~] = size(DATA);  % number of data sets per state and state variables
n_X     = size(X0mat,1)-1; % take number of states from X0mat instead of data set
ndata   = n_D * n_X;

if exist('DATAx','var') && ~isempty(DATAx)
    n_X2 = size(DATAx,2); % extra uni-variate dataset
else
    n_X2 = 0;
end

add_data=0; % check if empty data sets have been added
for i = 1:n_X % catch cases where data were forgotten ...
    if size(DATA,2) < i
        DATA{1,i} = 0;
        add_data = 1;
    end
end
for i = 1:ndata % run through state variables
    if isempty(DATA{i}) % catch cases where data were forgotten ...
        DATA{i} = 0;
        add_data = 1;
    end
end
if n_D > 1 && add_data == 1 % data have been added, so a reshape is needed
    DATA = reshape(DATA,n_D,n_X);
end
[n_D,n_X] = size(DATA);  % number of data sets per state and state variables
ndata     = n_D * n_X; % recalculate ndata
    
glo2.n_X  = n_X;
glo2.n_X2 = n_X2;
glo2.n_D  = n_D;

% if axis and legend labels are forgotten, just give default names
if ~isfield(glo,'ylab') || isempty(glo.ylab{1})
    ylab{1} = 'State var. 1';
end
if n_X ~= size(glo.ylab,2)
    for i = 2:n_X
        if length(glo.ylab) < i
            glo.ylab{i} = ['State var. ',num2str(i)];
        end
    end
end
if ~isfield(glo,'xlab')
    glo.xlab = 'time';
end
if ~isfield(glo,'leglab1')
    glo.leglab1 = 'Scen. '; 
end
if ~isfield(glo,'leglab2')
    glo.leglab2 = ''; 
end

% Special options for call_deri (not used in all packages!)
if ~isfield(glo,'useode')
    glo.useode = 1; % calculate model using ODE solver (1) or analytical solution (0)
end
if ~isfield(glo,'eventson')
    glo.eventson = 0; % events function for ODE solver on (1) or off (0)
end
if ~isfield(glo,'stiff')
    glo.stiff = 0; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
end
if ~isfield(glo,'break_time')
    glo.break_time = 1; % break time vector up for ODE solver (1) or don't (0)
end

if n_X ~= size(X0mat,1)-1 % this should not occur anymore
    error('The number of state variables in X0mat does not match the number of data sets entered. If you do not have data for a state, use DATA{x}=0.')
end
   
%% Preliminary checks on data set(s)

for i = 1:ndata % run through state variables
    if size(W,2) >= i
        sD = size(DATA{i}); % calculate size of data set i
        sW = size(W{i});    % calculate size of weight factors for set i
        if sD(1) == sW(1) && sD(2) == sW(2) % if they are exact same size ...
            W{i} = W{i}(2:end,2:end); % assume that time and scenario where added to weights, so trim them off
            sW = size(W{i});   % calculate size of weight factors for set i again
        end
        if sum(sW) ~= 0 && (sD(1)-1 ~= sW(1) || sD(2)-1 ~= sW(2)) % if they NOT are exact same size ...
            error(['in data set ',num2str(i),' the weight matrix does not match the data set.'])
        end
    end
    if DATA{i}(1,1) == -1 % if we have survival data ...
        a = DATA{i}(2:end,2:end); % extract data set
        b = diff(a); % deaths in each interval
        if any(b(:)>0) % make sure there are no zombies in there ...
            error(['in data set ',num2str(i),' number of survivors should never increase in time.'])
            % warning(['in data set ',num2str(i),' number of survivors should never increase in time.'])
            % % you can replace this error with a warning if you know what you're doing
        end
        if size(W,2) >= i
            w = W{i}; % extract weights
            if any(w(isnan(a))) > 0
                error(['for state number ',num2str(i),' do not fill in missing/removed animals at a point where the observation is NaN.'])
            end
        end
        mat_hlp = ones(n_D,1) * (1:n_X); % helper matrix to find which state we deal with
        X_hlp = mat_hlp(i); % this is the state we deal with (needed when there are multiple data sets!)
        if any(X0mat(X_hlp+1,:) > 1) % and check that the initial values are reasonable
            error(['for state number ',num2str(i),' initial values in X0mat should be probabilities (so <1).'])
        end
    end
    if DATA{i}(1,1) == -3 % if we have immobility data ...
        if size(W,2) >= i && any(W{i}>0)
            error(['for state number ',num2str(i),' I cannot deal (yet) with missing/removed animals for immobility analysis.'])
        end
    end
end

if isfield(glo,'wts') % check if weights are right size
    if sum(size(DATA) == size(glo.wts)) ~= 2 % sizes do not match
        error(['Weights matrix glo.wts must be same size as data sets matrix DATA: ',num2str(n_D),'x',num2str(n_X)])
    end
end
if isfield(glo,'var') % check if variances are right size
    if sum(size(DATA) == size(glo.var)) ~= 2 % sizes do not match
        error(['Variance matrix glo.var must be same size as data sets matrix DATA: ',num2str(n_D),'x',num2str(n_X)])
    end
end

if n_D > 1 % some extra care needed when there are additional data sets
    if numel(W) < numel(DATA) % W matrix is smaller than DATA
        W_new = cell(n_D,n_X); % create empty cell array with same size as DATA
        W_new(1:size(W,1),1:size(W,2)) = W; % enter the W array in there
        W = W_new; % replace old array with new one
    end
end

% Make sure that weight matrix is specified correctly
add_wts = 0; % keep track whether weights are added
for i = 1:ndata % loop over all data sets
    if isempty(W) || length(W(:))<i || isempty(W{i}) || numel(W{i}) < 2 % if weights are not specified for this data set ...
        if DATA{i}(1,1) == -1 || DATA{i}(1,1) == -3 % if we have survival data ...
            W{i} = zeros(size(DATA{i})-1); % make zeros with the same size as the data (without times and scenarios)
            % Note that for survival, the weights matrix specifies missing/removed animals
            add_wts = 1;
        else
            W{i} = ones(size(DATA{i})-1); % make ones with the same size as the data (without times and scenarios)
            add_wts = 1;
        end
        if sum(size(W{i}) == size(DATA{i})-1) ~= 2 % size of weight matrix does not match data set
            error(['Weight matrix of data set ',num2str(i),' is not the size of the data set.'])
        end
    end
end
if n_D > 1 && add_wts == 1 % weights have been added, so a reshape might be needed
    W = reshape(W,n_D,n_X);
end

% an initial survival probability different from one is meaningless and
% will lead to problems in the analysis (in the GUTS package, locS is used
% to specify the position of survival probability in the state vector)
if isfield(glo,'locS') && any(X0mat(1+glo.locS,:) ~= 1)
    error('Initial survival probablity needs to be one (it is not one in X0mat)!')
end

%% Preliminary checks on extra data set(s)

% we also need to check DATAx and Wx, but for now, fix Wx as done for W
if n_X2 > 0
    for i = 1:n_X2 % run through extra data sets
        if isempty(DATAx{i}) % catch cases where data were forgotten ...
            DATAx{i} = 0;
        end
        if size(Wx,2) >= i
            sD = size(DATAx{i}); % calculate size of data set i
            sW = size(Wx{i});   % calculate size of weight factors for set i
            if sD(1) == sW(1) && sD(2) == sW(2) % if they are exact same size ...
                Wx{i} = Wx{i}(2:end,2:end); % assume that time and scenario where added to weights, so trim them off
                sW = size(Wx{i});   % calculate size of weight factors for set i again
            end
            if sum(sW) ~= 0 && (sD(1)-1 ~= sW(1) || sD(2)-1 ~= sW(2)) % if they NOT are exact same size ...
                error(['in extra data set ',num2str(i),' the weight matrix does not match the data set.'])
            end
        end
    end
    % Make sure that weight matrix is specified correctly
    add_wts = 0; % keep track whether weights are added
    for i = 1:n_X2 % loop over all extra data sets
        if isempty(Wx) || length(Wx(:))<i || isempty(Wx{i}) || numel(Wx{i}) < 2 % if weights are not specified for this data set ...
            Wx{i} = ones(size(DATAx{i})-1); % make ones with the same size as the data (without times and scenarios)
            add_wts = 1;
            if sum(size(Wx{i}) == size(DATAx{i})-1) ~= 2 % size of weight matrix does not match data set
                error(['Weight matrix of extra data set ',num2str(i),' is not the size of the data set.'])
            end
        end
    end
end

% if axis and legend labels are forgotten, just give default names
if ~isfield(glo,'ylab2') || isempty(glo.ylab2{1})
    glo.ylab2{1} = 'Y-values extra data 1';
end
if n_X2 ~= size(glo.ylab2,2)
    for i = 2:n_X2
        if length(glo.ylab2) < i
            glo.ylab2{i} = ['Y-values extra data ',num2str(i)];
        end
    end
end
if ~isfield(glo,'xlab2') || isempty(glo.xlab2{1})
    glo.xlab2{1} = 'X-values extra data 1';
end
if n_X2 ~= size(glo.xlab2,2)
    for i = 2:n_X2
        if length(glo.xlab2) < i
            glo.xlab2{i} = ['X-values extra data ',num2str(i)];
        end
    end
end

%% Other checks on inputs, and create time/conc vectors

names      = fieldnames(par); % extract all field names of par (global)
ind_fittag = ~strcmp(names,'tag_fitted');
names      = names(ind_fittag); % make sure that the fit tag is not in names
glo2.names = names; % put it in the global
if sum(ind_fittag) ~= length(ind_fittag) % then there is a tag_fitted in par
    par = rmfield(par,'tag_fitted'); % and remove the tag (assume that calling this script is always BEFORE fitting)
end

% Create time and concentration vectors for fitting. This is now done here,
% always, because that allows the user to make a profile likelihood after a
% simulation.

if ~isfield(glo,'Tinit')
    glo.Tinit = 0; % if not defined, assume that t=0 is the start of the calculations
end
tvect = [glo.Tinit]; % initialise time vector from data set with initial time
cvect = []; % initialise concentration vector as empty matrix
for i = 1:ndata
    tvect = [tvect;DATA{i}(2:end,1)];  % append all time vectors from the data set
    % and make sure that the initial time is in there!
    cvect = [cvect;(DATA{i}(1,2:end))']; % append all concentration vectors from the data set
    lam = DATA{i}(1,1); % how to transform the data
    if isnan(lam) % if some smartass puts a NaN there ...
        error(['The first number in the data matrix for state ',num2str(i),' should be a number. Check the byom script for the meaning of that number.'])
    elseif (lam < 0 && ~ismember(lam,[-1 -2 -3])) || lam > 2
        error(['The first number in the data matrix for state ',num2str(i),' cannot be interpreted. Check the byom script for the meaning of that number.'])
    end
end
ttot = unique(tvect(tvect>=glo.Tinit)); % total vector from data set with unique time points (global)

if length(ttot) < 2 % if we have zero or only 1 time point, there is no data to use
    if ~isfield(glo,'t') % so there must be a glo.t specified
        error('If you do not enter data, make sure to define glo.t manually in your script under the heading: Time vector and labels for plots.')
    else % we can use glo.t instead
        ttot = glo.t; % 
    end
end

% some adaptations if the time vector of the data is really short! (needed
% for ODE solver, otherwise it returns too many time points)
if length(ttot) == 2
    ttot = [ttot(1) mean(ttot) ttot(2)];
end
glo2.ttot = ttot;

if ~isfield(glo,'t') % if it is not specified, base it on t in data
    glo.t = linspace(glo.Tinit,max(ttot)*1.02,100); % this vector is used for plotting
elseif glo.t(1) ~= glo.Tinit % if the long time vector does not match the initial time ...
    glo.t = linspace(glo.Tinit,max(glo.t),length(glo.t));
    % Create a new vector of same length that starts at glo.Tinit as we
    % have to start the model at the correct values for the initial states!
    % The time vector in the data set is checked in calc_start.
end    
glo.t = glo.t(:); % make sure it is a column vector (NEW)

if ~isfield(glo,'saveplt')
    glo.saveplt = 0; % by default, don't save any plots
end

ctot = X0mat(1,:);    % this is the concentration vector that will be used!
% Note: the model is run ONLY with the concentrations/scenarios that are
% given in the initial state matrix X0mat! Scenarios that are in the data
% set but not in X0mat are thus ignored.

if length(ctot) ~= length(unique(ctot))
    error('The scenario numbers in X0mat (first row) must be unique!')
end
glo2.ctot = ctot;

%% Preliminary checks on model parameters

all_ones = 0;
for i = 1:length(names) % run through all parameters
    a = eval(['par.',names{i}]);
    if length(a) == 1
        all_ones = all_ones + 1; % count number of parameters with a 1-vector
    else
        if length(a) < 4
            error(['par.',names{i},': check parameter definition (specify a vector with 4-5 values for each parameter).'])
        end
        if length(a) == 4
            a(5) = 1; % add a one (normal scale)
            eval(['par.',names{i},'=a;']);
        end
        if a(5) == 0 % look at log-fitted parameters
            if a(1) == 0 % initial value is zero
                error(['par.',names{i},': you fit one parameter on log scale and have an initial value of zero!'])
            end
            if a(3) == 0 % lower bound is zero
                error(['par.',names{i},': you fit parameter on log scale and lower range is zero!'])
            elseif a(4) == 0 % upper bound is zero
                error(['par.',names{i},': you fit parameter on log scale and upper range is zero!'])
            end
        end
    end
end

% Wrap the glo and glo2 structures into a single wrapper (this is an
% adaptation for the use of the parallel toolbox)
WRAP.glo  = glo;
WRAP.glo2 = glo2;

pmat   = packunpack(1,par,0,WRAP);  % transform structure into a regular matrix
if all_ones == length(names) % user specified a 1-vector for the parameters
    par = packunpack(2,0,pmat,WRAP);  % transform regular matrix into a structure
    % this is handy, as it turns par to a standard setting as 5-column matrix
elseif all_ones>0 % some where 1 and others longer vectors
    error('Make sure the parameters all have the same size vector.')
end
% if any parameter is outside the bounds, error
if any(pmat(:,1)<pmat(:,3) | pmat(:,1)>pmat(:,4))
    error('Initial value of one or more parameters is out of min-max bounds!')
end

%% Final things

% Make sure that versions of Matlab 2017a and newer do not automatically
% update figure legends when new info is plotted in the same plot. This
% cannot be set for older versions (they don't update the legend anyway).
% I first modified the legend statements in the plotting routines to
% prevent updating, but this yields errors with older versions.
if ~verLessThan('matlab','9.2')
    % -- Code to run in MATLAB R2017a and later here --
    set(groot,'defaultLegendAutoUpdate','off');
end

% clear the extra variables that were defined in this script
clear a b n_X n_X2 n_D sD sW add_wts i ttot ctot cvect tvect