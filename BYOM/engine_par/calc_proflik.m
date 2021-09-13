function [par_best,XingS,prof_coll,loglik_disp] = calc_proflik(par,parname,opt_prof,varargin)

% Usage: [par_best,XingS,prof_coll,loglik_disp] = calc_proflik(par,parname,opt_prof,varargin)
%
% Calculates a profile likelihood for parameter with the name <parname>
% from your parameter structure. Note that <parname> can be a single
% string, or a cell array of strings to profile multiple parameters. Use
% the string 'all' to profile all fitted parameters.
% 
% <par> is the total parameter structure, resulting from the optimisation.
% Various options can be set in <opt_prof> (see <prelim_checks.m>), a.o. to
% specify the level of detail of the profiling, detailed (1) or coarse (2,
% default). Optionally, as extra input, <opt_optim> to allow for refitting
% when a better optimum is found during profiling.
% 
% Optional outputs: edges of the 95% CI in <XingS>, better parameter
% estimates in <par_best> (might be only very slightly better), entire
% profile in <prof_coll> (log pars on log scale), and better log-likelihood
% if a better value was found (<loglik_disp>).
%
% This function applies a variable stepsize algorithm. When the likelihood
% ratio changes very little (less than <Lcrit_min>), the stepsize is
% increased (up to a maximum, specified by <Fstep_max>). When the lik.
% ratio changes too much (more than <Lcrit_max>), the algorithm tries again
% with a smaller stepsize (also bound to a minimum: <Fstep_min>). Note that
% the stepsize is used as a fraction of the parameter value that is tried.
% To prevent very small stepsizes when the value goes towards zero (as can
% be the case for effect thresholds), I also apply an *absolute* minimum
% stepsize (<Fstep_abs>), which is specified as a fraction of the best
% parameter value (<Xhat>) (unless it is zero, then algoritm takes
% something small!).
%
% Author     : Tjalling Jager
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2
global DATA W DATAx Wx X0mat

% Create the wrapper to allow the parallel toolbox to function.
WRAP.DATA  = DATA;
WRAP.W     = W;
WRAP.DATAx = DATAx;
WRAP.Wx    = Wx;
WRAP.X0mat = X0mat;
WRAP.glo   = glo;
WRAP.glo2  = glo2;

verbose  = opt_prof.verbose; % set to 0 to suppress output to screen, 2 to only suppress plots
brkprof  = opt_prof.brkprof; % set to 1 to stop profiling when better optimum is located, 2 to re-fit
saved    = opt_prof.saved; % set to 1 to plot/display saved profiles

if ~isempty(varargin) && ~isempty(varargin{1})
    opt_optim = varargin{1}; 
elseif brkprof == 2
    brkprof = 1; % we cannot refit without opt_optim
    if verbose ~= 0
        warning('off','backtrace')
        warning('No re-fitting will be done as opt_optim is left empty as input. Will break instead')
        warning('on','backtrace')
    end
end

names    = glo2.names;
filenm   = glo.basenm;

% predefine outputs, which may not always be set
loglik_disp = [];

pmat = packunpack(1,par,0,WRAP);  % transform structure into a regular matrix
if size(pmat,2) < 4 % no full parameter definition with bounds
    error('Specify the parameters in your script file with upper and lower bounds.')
end
if size(pmat,2) == 4 % when not fitted first, the 5th column is not added if needed
    pmat(:,5) = 1; % by default, assume normal scale
end
par_sel = pmat(:,2);      % this is the selection vector
ind_fit = pmat(:,2) == 1; % logical indices to fitted parameters

if ~iscell(parname) % check the input parameter names
    if strcmp(parname,'all') % an input 'all' implies all fitted parameters
        profpars = names(ind_fit);
    else % otherwise, the user just wants a single parameter
        profpars = {parname};
    end
else % otherwise, the input is a cell array of names
    if length(parname) == 1 && strcmp(parname{1},'all') % an input 'all' implies all fitted parameters
        profpars = names(ind_fit); % this catches the case where someone enters {'all'}
    else
        profpars = parname;
    end
end
n_pr = length(profpars); % number of profiles to make

% derive the location of the parameters in pmat (and thus in names)
parnum = zeros(n_pr,1);
for i = 1:n_pr
    parnum_tmp = find(strcmp(names,profpars{i})); % lookup the name in the structure and return the position
    if isempty(parnum_tmp)
        error(['The name of the parameter to profile (',profpars{i},') does not exist in your structure.'])
    else
        parnum(i) = parnum_tmp;
    end
end

if verbose ~= 0 && ~isfield(par,'tag_fitted') % apparently, parameters have not been fitted
    warning('off','backtrace')
    warning('Parameters have not been fitted, so results may not be very meaningful!')
    disp(' '), warning('on','backtrace')
end

if saved == 1 % use saved set
    if exist([filenm,'_LP.mat'],'file') ~= 2
        error(['There is no likelihood-profile saved with name: ',filenm,'_LP.mat'])
    end
    parnum_tmp = parnum; % remember the parnum we derived as we'll replace it with the saved set
    load([filenm,'_LP'])
    if ~isequal(parnum,parnum_tmp)
        error('The saved set contains different set(s) of fitted parameters that the requested ones.')
    end
    pmat = packunpack(1,par,0,WRAP);  % take pmat from the saved par, rather than the input
    par_best = [];
end

%% Set cut-off criteria for profiles based on chi-square

% list of critical values of the chi-square for 95% with df 1 to 5; that
% way we can work without the statistics toolbox in most cases
chitable = [3.8415 5.9915 7.8147 9.4877 11.07];

% For finding the bounds of the joint hypercube region, use number of
% parameters as dfs. This is used for the likelihood-region shooting method.
chicritJ = chitable(min(5,sum(par_sel==1))); % cut-off for the 95% parameter region
chicritS = chitable(1); % cut-off for single par. CI and propagation region

%% Initialise figure to plot the profiling ...

if verbose == 1
    
    % Calculate size of multiplot
    n = ceil(sqrt(n_pr));
    m = ceil(n_pr/n);
    [figh,ft] = make_fig(m,n); % create a figure window of correct size
    
    % collect plotting handles for later use
    h1_rem = cell(n_pr,1);
    
    for i = 1:n_pr
        h1 = subplot(m,n,i); % subplot
        hold on
        set(h1,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        % title(['Called from: ',scriptnm, ' (',date,')'],'Interpreter','none',ft.name,ft.title)
        
        h1_rem{i} = h1; % remember subplot handle
        
        % if the parameter was on log-scale, mention it in the axis label
        if pmat(parnum(i),5) == 0
            xlabel(h1,['parameter ',names{parnum(i)},' (log-scale)'],ft.name,ft.label)
        else
            xlabel(h1,['parameter ',names{parnum(i)}],ft.name,ft.label)
        end
        ylabel(h1,'minus 2x log-likelihood ratio',ft.name,ft.label)
        
        drawnow
    end
    
end

%% Calculate the profile(s)
% This uses the parallel toolbox. The up and down branches for all
% parameters are done in parallel. 

if saved ~= 1 % only when not using saved set
    
    % Start/check parallel pool
    if glo2.n_cores > 0
        poolobj = gcp('nocreate'); % get info on current pool, but don't create one just yet
        if isempty(poolobj) % if there is no parallel pool ...
            parpool('local',glo2.n_cores) % create a local one with specified number of cores
        end
    end
    
    disp('Profiling is done using the parallel toolbox. No progress will be shown!')
    disp(' ')
    
    run_profs = allcomb(parnum,[1 2]); % all combinations of parameters and up-down
    
    % Now, there is a while loop. If a better optimum is found there is an
    % automatic re-fit, and the profiling is redone. This probably needs an
    % option ... but then we can remove auto_profiles as a function.
    
    hurrah   = 0; % this will go to 1 if we can stop the while loop
    no_rnds  = 0; % count nr of times that we found a better optimum
    par_best = []; % this will catch the parameters for the best new optimum, if one is located
    
    while hurrah == 0
        
        % initialise cell arrays to capture profiling output
        Xcoll       = cell(size(run_profs,1),1);
        loglikrat   = Xcoll;
        par_better  = Xcoll;
        loglikmax   = Xcoll;
        sample_prof = Xcoll;
        
        parfor j = 1:size(run_profs,1) % calculate from profile down and up from the ML estimate
            [Xcoll_tmp,loglikrat_tmp,par_better_tmp,loglikmax_tmp,sample_acc_tmp] = ...
                sub_proflik(pmat,run_profs(j,1),run_profs(j,2),opt_prof,WRAP);
            Xcoll{j}      = Xcoll_tmp;
            loglikrat{j}  = min(loglikrat_tmp,1000); % this should not be needed, but sometimes we get an Inf in there ...
            par_better{j} = par_better_tmp;
            loglikmax{j}  = loglikmax_tmp;
            sample_prof{j} = sample_acc_tmp;
        end
        if no_rnds == 0 % in the first round ... collect the loglik
            loglikbest = -loglikmax{1};
        end
        
        % Next, find the overall best new optimum (if better optima are located),
        % and collect all 'acceptable' parameter sets for the likelihood-region
        % method.
        logliknew       = -0.01; % this ignores tiny, irrelevant, improvement
        sample_prof_acc = [];
        flag_better     = 0;
        for j = 1:size(run_profs,1)
            if min(loglikrat{j}) < logliknew
                logliknew = min(loglikrat{j}); % update new minimum MLL-ratio
                par_best  = par_better{j};     % update the better parameter set
                flag_better = 1;
                loglik_disp = loglikbest + logliknew/2; % this loglik_disp will only be used if we don't do extra optimisation!
            end
            sample_prof_acc = cat(1,sample_prof_acc,sample_prof{j}); % collect profile points in sample for likreg
        end
        
        % Next, a check whether we found a better optimum. A new fit is done,
        % but only if the profiled parameter is a fitted one! Otherwise, we'll
        % be stuck in a never-ending loop. This is a place for an option!
        if brkprof == 2 && flag_better == 1 && ind_fit(run_profs(j,1)) % the profiling has found a better value
            no_rnds = no_rnds + 1; % count a better optimum
            % update the new best set with a new optimisation
            opt_optim.it = 0; % don't show iterations of the simplex optimisation
            if verbose ~= 0
                warning('off','backtrace')
                warning('Better optimum was found, initiating new optimisation')
                warning('on','backtrace')
            end
            % Note that calc_optim will still put output to screen ...
            [par_best,loglik_disp] = calc_optim(par_best,opt_optim); % new par_tmp
            pmat     = packunpack(1,par_best,0,WRAP);  % transform structure into a regular matrix
        else
            hurrah = 1; % we can/must stop!
        end
        
    end
end

%% Check if we need to stop here already

if ~isempty(par_best) && brkprof == 1
    % if the profiler found a better optimum, and we need to break, break and return!
    diary (glo.diary) % collect screen output in the diary "results.out"
    disp(' ')
    fprintf('The profiles found a better optimum: %1g (best was %1g). \n',loglik_disp,loglikbest);
    disp('The profiles are not complete (the analysis is set to stop at a better optimum)')
    disp('   It is a good idea to also make new plots for the fit.')
    fprintf('Parameter values for new optimum (formatted for copy-paste) \n');
    fprintf('=================================================================================\n');
    print_par(par_best); % write the better estimate to screen in a formatted way
    fprintf('=================================================================================\n');
    diary off  % close results.out
    
    % initialise cell arrays to capture profiling output
    prof_coll  = cell(n_pr,1);
    Xhat       = nan(n_pr,1);
    
    % Since the parallel version does not plot on the fly, there is no
    % plotting in the graph. We can plot here what has been collected
    % before the break set in.
    for i = 1:n_pr % run through the profiles
        % find the correct indices in run_profs
        ind1 = run_profs(:,1)==parnum(i) & run_profs(:,2)==1;
        ind2 = run_profs(:,1)==parnum(i) & run_profs(:,2)==2;
        prof_coll{i} = [Xcoll{ind1}' loglikrat{ind1}';Xcoll{ind2}' loglikrat{ind2}']; % construct entire profile
        plot(h1_rem{i},prof_coll{i}(:,1),prof_coll{i}(:,2),'k.')
        Xhat(i) = pmat(parnum(i),1); % the best estimate for the parameter (2nd column is select)
        if pmat(parnum(i),5) == 0 % then we plot on log-scale
            Xhat(i) = log10(Xhat(i));
        end
        plot(h1_rem{i},Xhat(i),0,'ko','MarkerFaceColor','y') % plot the max lik value as a circle
    end
    
    XingS = [];
    return
end

%% Plotting the profiles and calculating the 95% intervals

% start by reconstructing prof_coll
if saved ~= 1 % only when not using saved set
    
    % initialise cell arrays to capture profiling output
    prof_coll  = cell(n_pr,1);
    
    for i = 1:n_pr % run through the profiles
        
        % find the correct indices in run_profs
        ind1 = run_profs(:,1)==parnum(i) & run_profs(:,2)==1;
        ind2 = run_profs(:,1)==parnum(i) & run_profs(:,2)==2;
        
        prof_coll{i} = [Xcoll{ind1}' loglikrat{ind1}';Xcoll{ind2}' loglikrat{ind2}']; % construct entire profile
        prof_coll{i} = sortrows(prof_coll{i},1); % sort on basis of x-axis
        % prof_coll will also be returned as output (log pars are still on log-scale)
    end
end

% next, calculate crossings, best fit set, and boundaries for calc_likregion to use
XingS      = cell(n_pr,1);
XingJ      = cell(n_pr,1);
boundscoll = zeros(n_pr,2); % initialise the matrix to catch the bounds of the hyperbox (for likreg)
Xhat       = nan(n_pr,1);

for i = 1:n_pr % run through the profiles
    
    Xhat(i) = pmat(parnum(i),1); % the best estimate for the parameter (2nd column is select)
    if pmat(parnum(i),5) == 0 % then we plot on log-scale
        Xhat(i) = log10(Xhat(i));
    end
    
    % The function calc_xing find the crossings with the critical values,
    % and also warns us if they run into the edges of the range.
    XingS{i} = calc_xing(prof_coll{i},chicritS);  % find crossings of the chicrit
    XingJ{i} = calc_xing(prof_coll{i},chicritJ);  % for total joint CR (or limited)
    boundscoll(i,:) = [XingJ{i}([1 end],1)]; % collect the extreme bounds for the confidence region
    % note that for the log-scale pars, boundscoll will stay on log scale!
    
    % back-transform Xing if needed for CIs on screen
    if pmat(parnum(i),5) == 0
        XingS{i}(:,1) = 10.^XingS{i}(:,1);
        XingJ{i}(:,1) = 10.^XingJ{i}(:,1);
    end
    
end

% make the bounds a bit larger, just to be on the safe side
boundscoll(:,1) = boundscoll(:,1)-0.05*(boundscoll(:,2)-boundscoll(:,1));
boundscoll(:,2) = boundscoll(:,2)+0.05*(boundscoll(:,2)-boundscoll(:,1));

if verbose == 1
    
    % plot the profiles into the already prepared figure window
    for i = 1:n_pr % run through the profiles
        
        % Plot the chi-square criterion for 95% and 1 degree of freedom
        minX = prof_coll{i}(1,1);
        maxX = prof_coll{i}(end,1);
        plot(h1_rem{i},[minX maxX],[chicritS chicritS],'k:')
        xlim(h1_rem{i},[minX-0.05*(maxX-minX) maxX+0.05*(maxX-minX)]) % limit x-axis
        
        % And make a nice plot of the profile
        plot(h1_rem{i},prof_coll{i}(:,1),prof_coll{i}(:,2),'k.-')
        
        % set the y-axis to allow to read the graph when lik-ratio becomes very large
        plotmin = min(prof_coll{i}(:,2));
        plotmax = max(10,min(30,max(prof_coll{i}(:,2)))); % force plotmax between 10 and 30
        ylim(h1_rem{i},[plotmin plotmax])
        
        plot(h1_rem{i},Xhat(i),0,'ko','MarkerFaceColor','y') % plot the max lik value as a circle

    end
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['profiles_',filenm];%
        save_plot(figh,savenm);
    end
end

%% Display on screen

if verbose ~= 0
    
    diary (glo.diary) % collect screen output in the diary "results.out"
    
    if ~isempty(par_best) % if the profiler found a better optimum, display it!
        disp(' ')
        fprintf('The profiles found a better optimum: %1g (best was %1g). \n',loglik_disp,loglikbest);
        if brkprof == 2
            disp('   Re-fitting was done, so the final results are for the better optimum.')
        end
        disp('   It is a good idea to also make new plots for the fit.')
        fprintf('Parameter values for new optimum (formatted for copy-paste) \n');
        fprintf('=================================================================================\n');
        print_par(par_best); % write the better estimate to screen in a formatted way
        fprintf('=================================================================================\n');
        
    end
    
    disp('  ')
    disp('95% confidence interval(s) from the profile(s)')
    disp('=================================================================================')
    
    for i = 1:n_pr % run through profiles
        
        if size(XingS{i},1) == 2 % no problems, it is a single interval
            if XingS{i}(1,2) == 1 % lowest crossing an upper boundary ...
                a = sprintf('<%#0.4g',XingS{i}(1,1));
            else
                a = sprintf(' %#0.4g',XingS{i}(1,1));
            end
            if XingS{i}(end,2) == 1 % highest crossing a lower boundary ...
                b = sprintf('>%#0.4g',XingS{i}(2,1));
            else
                b = sprintf(' %#0.4g',XingS{i}(2,1));
            end
            fprintf('%-10s interval: %10s - %s \n',names{parnum(i)},a,b)
        else
            fprintf('%-10s confidence interval is a broken set (check profile):\n',names{parnum(i)})
            for ix = 1:size(XingS{i},1)/2
                if ix == 1 && XingS{i}(1,2) == 1 % lowest crossing an upper boundary ...
                    a = sprintf('<%#0.4g',XingS{i}((ix-1)*2+1,1));
                else
                    a = sprintf(' %#0.4g',XingS{i}((ix-1)*2+1,1));
                end
                if ix == size(XingS{i},1)/2 && XingS{i}(end,2) == 1 % highest crossing a lower boundary ...
                    b = sprintf('>%#0.4g',XingS{i}((ix-1)*2+2,1));
                else
                    b = sprintf(' %#0.4g',XingS{i}((ix-1)*2+2,1));
                end
                fprintf('      interval %u: %9s - %s \n',ix,a,b)
            end
        end
        
    end
    disp('=================================================================================')
    
    disp(' ')
    disp(['Time required: ' secs2hms(toc)])
    diary off  % close results.out
    
end

% Save the profile information after profiling is finished. This allows us
% to use these results later for a likelihood-region analysis.
if ~isempty(par_best) % if the profiler found a better optimum, display it!
    par = par_best; % need to save the best
end
save([filenm,'_LP'],'par','par_sel','boundscoll','prof_coll','sample_prof_acc','parnum')
% sample_prof_acc is added as these will also be included into the random
% sample. The parnum is also saved. This may be used in the future to be
% certain about which parameters where actually profiled (it is possible to
% profile parameters that have not been fitted).


% =========================================================================
% =========================================================================


function [Xcoll,loglikrat,par_better,loglikmax,sample_prof_acc] = sub_proflik(pmat,parnum,profdir,opt_prof,WRAP)

% This sub-function does the actual profiling.

% read options from structure
prof_detail = opt_prof.detail; % detailed (1) or a coarse (2) calculation
if  ~ismember(prof_detail,[1 2])
    prof_detail = 2; % if a wrong value is entered, 'coarse' is default
end
n_sub         = opt_prof.subopt;  % number of sub-optimisations to perform to increase robustness
brkprof       = opt_prof.brkprof; % set to 1 to break the profiling when a better optimum is located
subann        = opt_prof.subann;  % set to 1 to use simulated annealing (followed by simplex) instead of suboptimisations
Fsub_opt      = opt_prof.subrng; % maximum factor on parameters (both higher and lower) for sub-optimisations

sub_opt     = [Fsub_opt-1/Fsub_opt 1/Fsub_opt];% [2.67 0.33]; % settings for the suboptimisation: par * (rand * sub_opt(1)+ subopt(2))
sub_opt_log = [(log10(sum(sub_opt))-log10(sub_opt(2))) log10(sub_opt(2))]; % modify to get same range for log parameters!

%% Initialise things

% put parameters that need to be on log scale on log10 scale
pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));
% also modify the min max ranges!
pmat2 = pmat; % but do this on a copy, as transfer wants the non-log version of the ranges
pmat2(pmat2(:,5)==0,[3 4]) = log10(pmat2(pmat2(:,5)==0,[3 4]));

parshat = pmat(:,1);       % this is the parameter vector
par_sel = pmat(:,2);       % this is the selection vector
ind_fit = par_sel == 1;    % logical indices to fitted parameters
Xhat    = pmat(parnum,1);  % the best estimate for the profiled parameter

par_better  = []; % start with empty output for this one; it will get filled when a better optimum is located

%% Set performance parameters
% These values can be changed if the performance of the likelihood
% procedure is not adequate.

% list of critical values of the chi-square for 95% with df 1 to 5; that
% way we can work without the statistics toolbox in most cases
chitable = [3.8415 5.9915 7.8147 9.4877 11.07];

% For finding the bounds of the joint hypercube region, use number of
% parameters as dfs. This is used for the likelihood-region shooting method.
chicritJ = chitable(min(5,sum(par_sel==1))); % cut-off for the 95% parameter region
% chicritS = chitable(1); % cut-off for single par. CI and propagation region

switch prof_detail
    case 1 % detailed options
        Fstep_min = 0.001;    % min stepsize (as fraction of value that is tried)
        Fstep_max = 30;       % max stepsize (as fraction of value that is tried)
        if Xhat == 0          % if parameter value is zero; than parameter is likely a NEC
            Fstep_abs = 1e-4; % This is just a low value that should be ok in most situations
        else
            Fstep_abs = 1e-3 * abs(Xhat); % smallest stepsize (in absolute sense)
        end
        % (this is to prevent problems when the parameter value that is tried goes
        % to zero; the stepsize is a fraction of the value tried so can also become
        % very small).
        
        % criteria for change in likelihood ratio
        Lcrit_max  = 1;   % maximum change in likelihood, above this value, stepsize is decreased
        Lcrit_min  = 0.2; % minimum change in likelihood, below this value, stepsize is increased
        Lcrit_stop = chicritJ+5;  % stop when likelihood ratio reaches this value
        
    case 2 % coarse options
        Fstep_min = 0.005;    % min stepsize (as fraction of value that is tried)
        Fstep_max = 30;       % max stepsize (as fraction of value that is tried)
        if Xhat == 0          % if parameter value is zero; than parameter is likely a NEC
            Fstep_abs = 1e-3; % This is just a low value that should be ok in most situations
        else
            Fstep_abs = 1e-2 * abs(Xhat); % smallest stepsize (in absolute sense)
        end
        % (this is to prevent problems when the parameter value that is tried goes
        % to zero; the stepsize is a fraction of the value tried so can also become
        % very small).
        
        % criteria for change in likelihood ratio
        Lcrit_max  = 2; % maximum change in likelihood, above this value, stepsize is decreased
        Lcrit_min  = 0.7; % minimum change in likelihood, below this value, stepsize is increased
        Lcrit_stop = chicritJ+2;  % stop when likelihood ratio reaches this value
        
end

% Set options for the fitting (rougher than for finding the best fit;
% increases speed). Also options are set to turn the iterations on screen
% off.
oldopts  = optimset('fminsearch');
options1 = optimset(oldopts,'Display','off','FunValCheck','on','TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',30*(length(parshat)-1)); % only do rough optimisation
options2 = optimset(oldopts,'Display','off','FunValCheck','on','TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50*(length(parshat)-1)); % do better optimisation

%% Start with the profiling ...

par_sel(parnum) = 0; % do not fit this parameter anymore, as we make profile
pmat_orig = pmat;    % remember the original parameter values (we are messing with a copy)
pmat(:,2) = par_sel; % overwrite the selection part of pmat with the new par_sel

p_index   = (par_sel==1); % find the parameters that need to be fitted in profiling
loglikmax = -1 * transfer(parshat(p_index),pmat,WRAP); % use transfer to obtain max likelihood!

parmin  = pmat2(parnum,3); % use min and max of parameter as provided in script
parmax  = pmat2(parnum,4); % taken from pmat2, so these may be on log-scale

flag_better  = 0; % keep track of whether a better value is found, so we can stop the analysis if needed
logliklowest = 0; % keep track of lower optima

switch profdir
    case 1 % calculate from profile down from the ML estimate
         proffact     = -1; % going down ...
         % plot(h1,Xhat,0,'ko','MarkerFaceColor','w') % plot the max lik value as a circle
    case 2 % calculate from profile up from the ML estimate
        proffact = 1; % going up ...
end

p_try     = pmat; % temporary parameter vector
Xtry      = Xhat; % initialise the parameter to try on the ML estimate
logliktmp = 0;
i         = 1;

% first values are zero
Xcoll(1)     = Xtry; % all tried parameter values will be stored in Xcoll
loglikrat(1) = 0;
flag         = 0;

% initialise Fstep as fraction of parameter value
Fstep   = Fstep_min; % start from the minimum step size
flagmin = 0; % test: flag is used to indicate that stepsize cannot or should not be decreased

sample_prof_acc = []; % initialise a matrix that will hold the points from the profile that will end up in the sample
% NOTE: the profiled, optimised, parameter sets will be added to the sample
% for the likelihood-region shooting method to increase robustness.

while flag == 0 % as long as we have not reached a stopping criterion
    
    flag_try = 0;
    % remember the values from the previous step in the profiling
    Xtry_old      = Xtry;
    logliktmp_old = logliktmp;
    
    while flag_try == 0 && flag == 0 % as long as we haven't got a good new value
        
        if Xtry_old ==0 % if old tryvalue is zero, we cannot use the relative stepsize.
            Xtry    = Xtry_old + proffact * Fstep_abs; % try a new value, based on absolute stepsize
            Fstep   = Fstep_min; % and re-initialise the relative stepsize
            flagmin = 1;
        elseif Fstep * abs(Xtry_old) >= Fstep_abs % new stepsize is larger than abs min size
            Xtry = Xtry_old + proffact * Fstep * abs(Xtry_old); % try a new value for the parameter
            % use the given relative stepsize, or the minimum absolute stepsize
            
            % now, if we are at least 1.5 times above the absolute
            % minimum and 1.5 times above the minimum relative stepsize
            % ... then it makes sense to allow for a decrease in
            % stepsize. Thus, flag is set to zero.
            if Fstep * abs(Xtry_old) > 1.5 * Fstep_abs && Fstep > 1.5 * Fstep_min
                flagmin = 0;
            else
                flagmin = 1;
            end
        elseif Fstep * abs(Xtry_old) < Fstep_abs
            Xtry    = Xtry_old + proffact * Fstep_abs; % try a new value for the parameter
            Fstep   = Fstep_abs / abs(Xtry_old);
            flagmin = 1;
        end
        
        % dont profile to negative values or very high ones!
        % make sure that values are within the specified parameter bounds
        % (bounds specified in the script file)
        Xtry = max(Xtry,parmin);
        Xtry = min(Xtry,parmax);
        
        p_try(parnum,1) = Xtry ; % put new profile value into the total try vector
        p_fit = p_try(p_index,1); % these are the ones that must be estimated now
        
        % Robustness of the optimisation routines is a real issue with
        % profiling, unfortunately. Here, I implement two methods: use
        % simulated annealling (probably less sensitive to local
        % minima), or a method with simplex sub-optimisations. The
        % latter seems to be best (but that will depend on the settings
        % of annealing).
        
        if subann == 1
            
            opt_ann.InitTemp  = 5; % much higher initial temperature
            opt_ann.Verbosity = 0; % controls output on screen (no output here)
            opt_ann.StopTemp  = 1e-2; % crude temperature for stopping as Simplex will finish it off
            [parshattemp,~] = anneal(@transfer,p_fit,opt_ann,p_try,WRAP); % first rough estimation
            
        else
            
            % Now do fminsearch to get new best fit, with profile parameter fixed
            % Do the optimisation at least twice: once rough, and once a little more
            % precise. Hopefully, this will avoid local minima, but still
            % be rapid enough!
            
            % minimisation is done with the reduced parameter set, other parameters in par are kept
            % to the fixed value. This means that parshat will also contain the limited set only.
            [parshatsub,logliksub] = fminsearch('transfer',p_fit,options1,p_try,WRAP);
            % Next, there is an option to perform additional
            % sub-optimisations with perturbed starting values. This
            % seriously increases robustness, but calculation time as well.
            indmodlog = p_try(p_index,5)==0; % index to log parameters for fitted parameters
            for i_sub = 1:n_sub
                modpar = p_fit .* (rand(length(p_fit),1)*sub_opt(1)+sub_opt(2));
                if any(indmodlog==1) % for log-parameters, we need a different approach to get the same range
                    modpar(indmodlog) = p_fit(indmodlog) + rand(sum(indmodlog),1)*sub_opt_log(1)+sub_opt_log(2);
                end
                % make sure they are all within their bounds, otherwise fminsearch cannot start
                modpar = max(modpar,pmat2(p_index,3));
                modpar = min(modpar,pmat2(p_index,4));
                
                [parshatTMP2,logliknew] = fminsearch('transfer',modpar,options1,p_try,WRAP);
                % modify starting point by randomly making parameters higher and lower
                logliksub(i_sub+1) = logliknew;
                parshatsub(:,i_sub+1) = parshatTMP2;
            end
            [~,locsub] = min(logliksub); % find lowest loglik value from the sub-sampling
            parshattemp = parshatsub(:,locsub); % and use the parameter values for that one
        end
        
        % use a Simplex to finish it
        [parshattemp,logliknew] = fminsearch('transfer',parshattemp,options2,p_try,WRAP);
        
        % and make the vector whole again
        p_try(p_index,1) = parshattemp ; % replace values in total vector with values that were fitted
        logliktmp        = -2 * (-1*logliknew-loglikmax); % likelihood ratio that follows a chi square
        deltaL           = abs(logliktmp - logliktmp_old); % difference with previous likelihood
        
        if logliktmp <= chicritJ % if the point is in the total 95% conf. region
            sample_prof_acc  = [sample_prof_acc ;[p_try(ind_fit,1)' logliktmp]]; % add it as 'accepted set'
        end

        % If the change in likelihood ratio is larger than we allow, or
        % if we are at the border of the allowable parameter range, we
        % want to decrease the stepsize. If the smalles stepsize is
        % already used, we can stop the analysis.
        skippit = 0;
        
        if deltaL > Lcrit_max
            if flagmin == 0 % & logliktmp < 5,
                % Only reject try and decrease stepsize when the flag says that it is possible!
                % Test: and only when we are in an interesting range of likelihood values!
                Fstep   = max(Fstep / 2 , Fstep_min); % try again with half the stepsize, unless less than min step
                skippit = 1; % no storage for this one!
            elseif any(Xtry == [parmin parmax])
                flag = 1; % if already using smallest stepsize, stop!
            end
        elseif any(Xtry == [parmin parmax]) % also stop when max or min is reached and deltaL < Lcrit_max!
            flag = 1; % stop!
        end
        
        % If the change is small enough, or if we want to stop the analysis anyway, store the value.
        if skippit == 0
            flag_try = 1;
            i = i+1;
            loglikrat(i) = logliktmp; % this is 2x log-likelihood ratio
            Xcoll(i)     = Xtry; % also remember the value of the profiles parameter
            
%             if verbose == 1
%                 % plot(h1,Xcoll(i),loglikrat(i),'k.') % plot em on the fly!
%                 drawnow
%             end
            
            if brkprof > 0 && flag_better == 1 && logliktmp > logliklowest
                % We have found a better value than the best one
                % already (flag_better=1), but now it's getting worse
                % again, so can break off the run.
                par_better = packunpack(2,0,pmat_better,WRAP); % return the better estimate in par_better
                return     % return to function calling us
            end
            
            if logliktmp < logliklowest % if we found a lower likelihood, remember it!
                if logliktmp < -0.01 % ignore tiny improvements before breaking the analysis
                    flag_better  = 1; % mark that we have located a better optimum
                end
                
                logliklowest = logliktmp;
                pmat_better  = [p_try(:,1) pmat_orig(:,2:5)]; % remember the better parameter set
                % put parameters that were on log scale back on normal scale
                pmat_better(pmat_better(:,5)==0,1) = 10.^(pmat_better(pmat_better(:,5)==0,1));
                par_better = packunpack(2,0,pmat_better,WRAP); % return the better estimate in par_better
                
            end
            
            if logliktmp > Lcrit_stop % stop if probability is high enough
                flag = 1;
            elseif deltaL < Lcrit_min % but change in likelihood is very small
                Fstep = min(Fstep * 2 , Fstep_max); % so double the stepsize, unless greater than max step
            end
        end
    end
end

