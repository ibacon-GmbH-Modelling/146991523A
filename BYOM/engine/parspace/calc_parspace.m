function [minmax,mll,stats,pmat_print,pmat] = calc_parspace(pmat,opt_optim,SETTINGS_OPTIM)

% Usage: [minmax,mll,stats,pmat_print,pmat] = calc_parspace(pmat,opt_optim,SETTINGS_OPTIM)
%
% Optimisation and finding a likelihood-based joint confidence region
% without the need for starting values. However, the parameter space to be
% searched must be limited for the algorithm to be effective. Parameter
% space to be searched is defined as the min-max bounds of the parameters
% in the matrix <pmat>. This is basically a genetic algorithm, whose aim is
% to find the best fitting parameter set and a cloud of points within a
% certain goodness-of-fit criterion (to be used for calculating confidence
% intervals). The algorithm applies likelihood profiling to refine and test
% the cloud from parameter space.
%
% Inputs
% <pmat>        parameter matrix
% <opt_optim>   structure with options for optimisations
% <SETTINGS_OPTIM> = settings read from setup_settings by <calc_optim_ps>
% 
% Outputs
% <minmax>  this is used to restart when slow kinetics is indicated; in
%    that case, it is a matrix with the min/max for each parameter of the
%    cloud that was used for testing slow kinetics. If all goes well,
%    <minmax>=-1 and no restart is made.
% <mll>     minus log-likelihood; this will be the lowest value found in the run
% <stats>   some statistics from the run: vector with number of sets in
%           joint CI, inner rim, and propagation band
% <pmat_print>  matrix with information on each parameter (best value and CI info)
% <pmat>        parameter matrix with best-fit values in first column
%
% Author     : Tjalling Jager
% Date       : June 2021
% Web support: <http://www.debtox.info/byom.html> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <calc_parspace> code that is
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

% Derive the <names> cell array since this will be saved with the MAT file
names = glo2.names; % also save the <names> variable with the sample
ind_fittag = ~strcmp(names,'tag_fitted');
names      = names(ind_fittag); % make sure that the fit tag is not in names_saved

%% BLOCK 1. Initial things and save when no parameters are fitted

% BLOCK 1.1. Make sure a few pieces of info are available if we want to
% save. We need to save very rapidly when there are no parameters to be
% fitted at all.

% debug_info    = 0; % when set to 1, prints extra information for debugging purposes
plot_intermed = opt_optim.ps_plots; % when set to 1, makes intermediate plots of parameter space to monitor progress
figh          = []; % start with an empty figure handle

% Use the fancy progress bar from Matlab
f = waitbar(0,'Initial things','Name','Calibration');

% BLOCK 1.2. % Derive several useful properties of <pmat>.
n_fit   = sum(pmat(:,2));       % number of fitted parameters
% n_par   = size(pmat,1);         % number of parameters in total
ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
minmax  = -1; % set to -1, so by default <calc_optim> thinks all is okay
ind_log_all = find(pmat(:,5) == 0); % indices to log-scale parameters in <pmat> (vector)
ind_log = find(pmat(ind_fit,5) == 0); % indices to FITTED log-scale parameters (vector)

% BLOCK 1.3. When no parameters are fitted, we need to do something. In
% principle, someone can fix all parameters and just run that (as
% simulation). That is fine, but since the plotting and post calculations
% all depart from a saved calibration file, we need to create one without
% going through the optimisation.
if isempty(ind_fit) % no fitted parameters at all!
    
    coll_all  = []; % make empty sample from parameter space
    coll_prof_pruned = []; % make empty profile likelihoods
    stats     = [0,0,0]; % fill the statistics (for <calc_optim> to print) with zeros
    
    pmat_print      = nan(size(pmat,1),8); % define a <pmat_print> with NaNs
    pmat_print(:,1) = pmat(:,1); % copy the input parameter values to the first column of <pmat_print>
    
    % Save a calibration file as MAT file.
    save([glo.basenm,'_PS'],'pmat','coll_all','pmat_print','coll_prof_pruned','names')
    disp_pmat(pmat,pmat_print); % display pmat in a formatted manner on screen
    
    mll             = transfer([],pmat); % calculate the min-log-likelihood for this parameter combination
    % Note: the MLL is not saved in the MAT file, but it is printed on
    % screen. When using a saved sample, the MLL cannot be reported and is
    % given as INF. I don't think that is an issue at the moment.
    
    delete(f) % delete progress bar again
    return % we can leave this function now and return to <calc_optim>
    
end

% BLOCK 1.4. Put <pmat> on log-scale where needed, and create a new matrix
% with bounds that are log-transformed where needed. Note: first column of
% <pmat> will be log-transformed where needed, which is already done in
% <calc_optim>. 
bnds_tmp = pmat(:,[3 4]); % copy parameter bounds to a temporary version
bnds_tmp(ind_log_all,[1 2]) = log10(bnds_tmp(ind_log_all,[1 2])); % and log10-transform these bounds where needed
bnds_tmp = bnds_tmp(ind_fit,:); % only keep the ones for the fitted parameters
% We want to keep the bounds in <pmat> itself on the normal scale. Reason
% is mainly that we'll save <pmat> at the end of this function and want to
% avoid rounding errors in the ranges from transforming back and forth.

% extract some settings for catching slow kinetics
if isfield(glo,'loc_kd') && isfield(glo,'loc_mi') % if these are defined, then we are doing a GUTS or DEBtox analysis
    % Note: location of parameters in total parameter vector is defined in
    % <startgrid_guts> and <startgrid_debtox>.
    loc_kd = glo.loc_kd;
    loc_mw = glo.loc_mi(1); % there can be more thresholds: use first one indicated in <startgrid_...> for checking slow kinetics
    % also find location of these parameters in the vector with FITTED parameters
    loc_kd_fit = find(ind_fit==loc_kd);
    loc_mw_fit = find(ind_fit==loc_mw);
    if pmat(loc_kd,2)==1 && pmat(loc_mw,2)==1 && pmat(loc_mw,5)==1 && pmat(loc_kd,5)==0 && bnds_tmp(loc_mw_fit,2)/bnds_tmp(loc_mw_fit,1) > 10
        % Only perform this check when both <kd> and <mw> are fitted, <mw>
        % is on normal scale, and <kd> is on log-scale. If one of these
        % things is different, the user has tampered with the parameter
        % ranges and he/she is on his/her own! And, also don't do this
        % check when the range of <mw> is less than a factor of 10 (when
        % the user has modified the range); for a small range, normal scale
        % will generally suffice.
    else
        loc_kd = -1; % to tell the slow-kin-catcher that it does not need to do anything
    end
else
    loc_kd = -1; % to tell the slow-kin-catcher that it does not need to do anything
end

%% BLOCK 2. Settings for the parameter space algorithm
% Many of the settings are globals, defined in <initial_setup>.

% Options for parameter optimisation. Note that many of these settings are
% vectors: they change with each subsequent round of the analysis.
crit_add = SETTINGS_OPTIM.crit_add; % extra on top of chi2-criterion to select ok values for next round
n_tr     = SETTINGS_OPTIM.n_tr; % number of random parameter tries in round 2 and further, for each ok parameter set
f_d      = SETTINGS_OPTIM.f_d;  % maximum step for random mutations, as factor of initial grid spacing
n_ok     = SETTINGS_OPTIM.n_ok(n_fit); % number of values that will at least continue to next round
n_conf   = SETTINGS_OPTIM.n_conf_all(n_fit,:); % stop criterion: minimum number of values within the total joint CI and inner rim, based on <n_fit>
tries_1  = SETTINGS_OPTIM.tries*ones(1,n_fit); % vector with number of tries per parameter for the initial grid
% NOTE: NOW TAKE SAME NUMBER FOR ALL PARAMETERS. TRIES_1 BECOMES A VECTOR WITH LENGTH OF TOTAL NR OF
% FITTED PARS IN PMAT (THIS IS DIFFERENT FROM THE OPENGUTS IMPLEMENTATION).
n_max    = SETTINGS_OPTIM.n_max;  % maximum number of rounds for the algorithm (default 12)

% We work here with the log-likelihood itself, so the chi2 criterion needs to be divided by 2.
chicrit_joint  = 0.5 * SETTINGS_OPTIM.crit_table(n_fit,1); % criterion for joint 95% CI or parameters
chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1);     % criterion for single-parameter CIs
chicrit_rnd    = chicrit_joint + crit_add;      % vector with criteria for selecting new tries in each subsequent round
chicrit_max    = max(1.2 * chicrit_joint,chicrit_single+1.5); % criterion for, roughly, a 97.5% joint CI for pruning the sample ...
% Note: for 1 or 2 fitted parameters, we need a bit more than 1.2 times the
% joint criterion, to make sure the profiles look nice and resampling is
% not triggered unnecessarily. Therefore the <chicrit_single+1.5> as minimum.

% Make sure these are defined, as we might end this function prematurely
% when hitting slow kinetics.
stats      = -1; % some statistics of the run
pmat_print = -1; % this matrix collects best value and CIs

%% BLOCK 3. Round 1 is using a regular grid over parameter space
% The first round is special and different from the later rounds. Create
% vectors with values to try for each parameter as a regular-spaced range
% between the min and max bounds (for log-scale parameters, the range is
% thus log-linear).
%
% Here, a cell array is used (<p_try>) as the number of trials differs
% between parameters. This array is turned into a regular matrix with all
% permutations of parameter values.

% BLOCK 3.1. Initialisation.
n_rnd   = 1;             % counter for rounds of optimisation
p_try   = cell(1,n_fit); % initialise <p_try> as empty cell array (one cell for each parameter)
d_grid  = nan(1,n_fit);  % initialise the grid spacing vector with NaNs (will collect spacing for each parameter)

% BLOCK 3.2. Create parameter vectors for each fitted parameter, with a regular grid.
for i_p = 1:n_fit % run through fitted parameters
    p_try{i_p}  = linspace(bnds_tmp(i_p,1),bnds_tmp(i_p,2),tries_1(i_p)); % vector: first tries as linear range between min-max bounds
    d_grid(i_p) = (bnds_tmp(i_p,2)-bnds_tmp(i_p,1))/(tries_1(i_p)-1); % difference between parameter values for parameter <i_p> (grid spacing)
end

% BLOCK 3.3. Create a large matrix <coll_all> with all permutations of the
% parameter values in <p_try>. For more than 5 fitted parameters, however,
% a regular grid is impossibly slow ... therefore use Latin-Hypercube
% sampling instead! If the user does not have the statistics toolbox,
% regular random sampling will be used, but this will be less effective
% (not tested). LHS sampling now also used for 5 parameters and rough
% settings.

if n_fit < 5 || (n_fit == 5 && opt_optim.ps_rough == 0)% use the default approach of a regular grid
    
    % This is the same as in openGUTS
    n_tries  = prod(tries_1); % total number of tries in this round (all permutations)
    coll_all = 1./zeros(n_tries,n_fit+1); % initialise matrix to catch all tries and their minloglik with INF
    
    coll_all(:,1:n_fit) = allcomb(p_try{1:n_fit}); % use smart function to make all permutations for all parameters (i.e., create a grid)
    % THIS ALSO WORKS WHEN NR OF TRIES DIFFERS FOR THE PARAMETERS!
    
else % THIS IS NEW AFTER V5.1 of BYOM
    
    n_tries  = 10000; % number of elements in the latin hypercube sample
    n_ok     = min(n_ok,500);   % test: default for 6/7 pars is 800, but for 5 it is 400
    coll_all = 1./zeros(n_tries,n_fit+1); % initialise matrix to catch all tries and their minloglik with INF
    
    if exist('lhsdesign','file')~=2 % when lhsdesign does not exist exists as an m-file in the path
        sample_lhs = rand(n_tries,n_fit); % uniform random sample between 0 and 1
        disp('Using random sampling in first round (Latin hypercube requires stats toolbox), instead of regular grid.')
    else
        sample_lhs = lhsdesign(n_tries,n_fit); % Latin-hypercube sample between 0 and 1
        disp('Using latin-hypercube sampling in first round, instead of regular grid.')
    end
    
    if n_fit > 5 && chicrit_joint < 6
        disp('Using a limited outer rim (based of df=5 rather than the actual number of parameters)')
        disp('   This implies that the blue points can no longer be interpreted as the joint CI!')
    end
    
    for i_p = 1:n_fit % go through the fitted parameters
        sample_lhs(:,i_p) = sample_lhs(:,i_p)*(bnds_tmp(i_p,2) - bnds_tmp(i_p,1))+bnds_tmp(i_p,1);
        % and change the sample to cover the bounds of the hypercube
    end
    coll_all(:,1:n_fit) = sample_lhs; % place sample in <coll_all>
    clear sample_lhs; % clear the sample to save memory
    
end

% BLOCK 3.4. Run through all elements of <coll_all> and calculate their
% likelihood. This could be integrated into the previous series of <for>
% loops ...
disp(' ')
disp(['Starting round 1 with initial grid of ',num2str(n_tries),' parameter sets'])
for i_t = 1:n_tries % run through all elements of <coll_all> and collect their likelihood
    waitbar(i_t/n_tries,f,'Round 1: initial grid') % update progress bar
    pfit              = coll_all(i_t,1:n_fit); % next, try this set of parameter values
    coll_all(i_t,end) = transfer(pfit,pmat); % calculate the min-log-likelihood for this parameter combination,
    % and collect it in the last column of <coll_all>
end

% BLOCK 3.5. Extract some useful matrices from the total <coll_all> matrix.
% Decide which sets to continue with in the next round. These will go into
% the new matrix <coll_ok>.
coll_all  = coll_all(~isinf(coll_all(:,end)),:); % remove the ones that have INF as minloglik
coll_all  = sortrows(coll_all,n_fit+1); % sort based on the minloglik in the last column (keep parameter sets together)
mll       = coll_all(1,end);            % lowest MLL; <coll_all> is sorted, so first is the best fitting one so far
ind_cont  = find(coll_all(:,end) < mll + chicrit_rnd(1),1,'last'); % index to last MLL that is within the criterium to continue with
ind_cont  = max(ind_cont,n_ok); % take at least the <n_ok> best ones ...
ind_cont  = min(ind_cont,size(coll_all,1)); % <n_ok> should always be smaller than the length of <coll_all>
% This latter check is only relevant when we fit two (or one?) parameters,
% otherwise the number of initial tries will always exceed <n_ok>.

% Check if we found a lot of ok values (that can for example happen when
% the ranges are set much tighter by the user).
if ind_cont > 1 * n_conf(1) % if we already have more than what we finally need ...
    ind_cont2 = find(coll_all(:,end)-mll > chicrit_rnd(2),1,'first'); % how many within *next* chi2 criterion?
    ind_cont  = max(n_conf(1),ind_cont2); % take highest from end number or the ones within the next chi-square criterion
    % This ensures that rather bad values (within <chicrit_rndi(2)>) are
    % still included at this point. We should take care not to remove too
    % many values-to-try early in the run.
end
coll_ok = coll_all(1:ind_cont,:); % take the <ind_cont> best ones to continue with

% BLOCK 3.6. Display status on screen, make plot, and prepare settings for next round.
disp(['  Status: best fit so far is (minloglik) ',num2str(mll)])

if plot_intermed == 1
    % And make a plot of the progress so far (one plot that will be updated after each round)
    figh = plot_grid(pmat,coll_ok,coll_all,[],figh,SETTINGS_OPTIM); % returns the handle to the graph in <figh>, so we can update the same plot
    uistack(f,'top')  % but place progress bar on top!
end

% Set all settings for the next round of optimisation.
n_rnd     = n_rnd + 1;   % increase counter for rounds by 1
n_tr_i    = n_tr(n_rnd); % number of random parameter tries in round 2
f_d_i     = f_d(n_rnd);  % maximum step as factor of grid spacing for random search
chicrit_i = chicrit_rnd(n_rnd); % chi2 criterion to select ok values

% Also check after first round if it's not too many. If we have a lot of
% parameter sets to continue with, we can use less tries in the next round.
if ind_cont > (1/2) * n_conf(1)   % if we have more than half of what we finally need ...
    n_tr_i = floor(n_tr_i/2);     % decrease the number of tries per set for next round
    if ind_cont >= 1 * n_conf(1)  % if we have more than what we finally need in total ...
        n_tr_i = floor(n_tr_i/2); % AGAIN decrease the number of tries per set
    end
end
n_tr_i = min(n_tr_i,max(2,floor(10*n_conf(1)/ind_cont))); % hard limit for the number of tries per set for next round
% This allows a max of 10x the target value to be tried in the next round
% (scaling back <n_tr_i>, with a minimum of 2). These checks should ensure
% that we don't have a huge amount of sets to try in round 2.

%% BLOCK 4ALT. Use slice sampling rather than the openGUTS mutations
% At this moment, this option is restricted to Tjalling. Since I modified a
% Matlab function, distributing it on the web is likely a copyright
% infringement.

flag_stop  = 0; % flag for when we can stop the analysis (sufficient points found: 1)

if opt_optim.ps_slice == 1 && exist('slicesample_byom','file')~=2 
    % only when slicesample_byom exists as an m-file in the path
    
    % We don't want to prune coll_all when going into the slice sampler;
    % pruning will be done there were needed.
    
    % Collect various input parameters for the slow-kinetics-catcher into a
    % structure.
    SLOKIN.loc_kd     = loc_kd;
    if loc_kd ~= -1
        SLOKIN.loc_mw_fit = loc_mw_fit;
        SLOKIN.loc_kd_fit = loc_kd_fit;
        SLOKIN.bnds_tmp   = bnds_tmp;
    end
    
    [coll_all,flag_stop,n_rnd,minmax] = calc_parspace_slice(coll_all,pmat,figh,plot_intermed,SLOKIN,SETTINGS_OPTIM);
    if flag_stop == 0 % then we returned prematurely, because slow kinetics was found
        return % so return to calc_optim_ps
    end
    % Note: under normal conditions, flag_stop=1 after the call to
    % calc_parspace_slice, which also implies that the regular mutation
    % rounds in BLOCK 4 will be skipped.
    
else % prepare for regular openGUTS mutation rounds
    
    % Prune <coll_all> to remove all values that are outside highest chi2 criterion.
    coll_all = coll_all(coll_all(:,end) < mll + chicrit_max,:);
    
end

%% BLOCK 4. Subsequent rounds are automated in a while loop

flag_inner = 0; % flag for when we will focus on inner rim (1)

while flag_stop ~= 1 % continue until this flag is set to 1
    
    disp(['Starting round ',num2str(n_rnd),', refining a selection of ',num2str(size(coll_ok,1)),' parameter sets, with ',num2str(n_tr_i),' tries each'])
    % Note: waitbar updates are now dealt with within <rand_mutations>

    % BLOCK 4.1. Call <rand_mutations> to mutate <coll_ok>, give each new set
    % an MLL, and return it in <coll_tries>.
    coll_tries = rand_mutations(pmat,coll_ok,bnds_tmp,n_tr_i,f_d_i*d_grid,f,n_rnd);
    coll_all   = cat(1,coll_all,coll_tries); % add the tries to the total <coll_all>
    coll_all   = sortrows(coll_all,n_fit+1); % sort the combined set based on minloglik
    
    % BLOCK 4.2. Do an optimisation here, with low precision, to improve
    % the best value so far. Add the optimised best set to <coll_all>.
    pfit = coll_all(1,1:n_fit)'; % copy best values to <pfit> (<coll_all> is sorted, so first is the best fitting one so far)
    % Only parameter values that are to be fitted; use as starting values for fitting.
    [phat,mll] = setup_simplex(pfit,0,pmat); % do a rough optimisation
    coll_all   = cat(1,coll_all(1,:),coll_all);  % copy the previous best to the first position
    coll_all(1,:) = [phat' mll]; % update the best fitting one with the new fitted parameters and MLL
    % Note: this is a bit awkward but ensures that also the fixed values
    % are correctly copied (<phat> only contains the fitted parameters).
    
    % BLOCK 4.3. Derive some useful indices from <coll_all> and
    % <coll_tries>. Find the last value in <coll_all> and <coll_tries> that
    % still is within a certain criterion.
    ind_final  = find(coll_all(:,end)   < mll + chicrit_joint,1,'last');  % index in <coll_all> for total joint CI
    ind_single = find(coll_all(:,end)   < mll + chicrit_single,1,'last'); % index in <coll_all> for inner rim
    ind_cont_t = find(coll_tries(:,end) < mll + chicrit_i,1,'last');      % index in just-tried sets that qualify for continuation to next round
    ind_cont_a = find(coll_all(:,end)   < mll + chicrit_i,1,'last');      % index in <coll_all> that qualify for continuation to next round
    ind_inner  = find(coll_all(:,end)   < mll + chicrit_single + 0.2,1,'last'); % index for inner rim, with a little extra
    
    % Two additional indices (roughly for a 97.5% joint CI) in case we run
    % into trouble (if we try too many new parameters with chicrit_i).
    % These are indices to the last entry in the two matrices that still
    % are within the total cloud that we like to calculate (slightly more
    % than <chicrit_joint>).
    ind_cont_t2 = find(coll_tries(:,end) < mll + chicrit_max,1,'last');
    ind_cont_a2 = find(coll_all(:,end)   < mll + chicrit_max,1,'last');
    
    disp(['  Status: ',num2str(ind_final),' sets within total CI and ',num2str(ind_single),' within inner. Best fit: ',num2str(coll_all(1,end))])

    % BLOCK 4.4. Try to catch slow kinetics at this point: if there are
    % signs, redo the first round with <mw> on log-scale! It tests on a
    % part of <coll_all>: only the sets that are within the joint 95% CI
    % (with a minimum of <n_ok> sets).
        
    if loc_kd ~= -1 % only when we identified that it needs to be checked
        
        coll_tst = coll_all(1:max(ind_final,n_ok),:); % parameter set used to test for slow kinetics
        
        min_mw = min(coll_tst(:,loc_mw_fit)); % minimum (non-zero) value of <mw> in the test matrix
        min_kd = min(coll_tst(:,loc_kd_fit)); % minimum value of <kd> in the test matrix
        
        check_corr = corrcoef(log10(coll_tst(:,loc_mw_fit)),coll_tst(:,loc_kd_fit)); % correlation between <mw> and <kd> on log scale
        % For slow kinetics, <mw> and <kd> will extend to the lower part of
        % their range. This might fail when users modify the ranges.
        % However, then we still have the check on the correlation
        % coefficient.
        crit_mw = (min_mw - bnds_tmp(loc_mw_fit,1)) / diff(bnds_tmp(loc_mw_fit,[1 2])); % distance from lower bound as fraction of range
        crit_kd = (min_kd - bnds_tmp(loc_kd_fit,1)) / diff(bnds_tmp(loc_kd_fit,[1 2])); % distance from lower bound as fraction of range
        
        % The check is now pretty strict ... too strict? Want to avoid
        % putting <mw> on log scale when we run into single-dose runaway ...
        if check_corr(2) > SETTINGS_OPTIM.slowkin_corr % is there a sufficiently strong positive correlation between <mw> and <kd>?
            if crit_mw < SETTINGS_OPTIM.slowkin_pars || crit_kd < SETTINGS_OPTIM.slowkin_pars
                % And is either <kd> or <mw> at lower part of their bounds?
                close(figh) % close the figure we have been building
                delete(f) % delete progress bar
                
                % Calculate min and max, of the parameter cloud we tested, for
                % each parameter (this creates a matrix). This will be used to
                % create new (smaller) starting ranges for a new round.
                minmax = [min(coll_tst(:,1:end-1),[],1)' max(coll_tst(:,1:end-1),[],1)'];
                minmax(ind_log,:) = 10.^minmax(ind_log,:); % and put on normal scale where needed
                return % and go back to <calc_optim_ps> to force a restart
            end
        end
        
    end
    
    % BLOCK 4.5. This is the tricky bit! This block is aiming to make sure
    % not to try too few or too many points in the next round, but select
    % an efficient set for continuation in <coll_ok>.
    %
    % Note: the settings here are all tweaked to get good results in the
    % cases tested. These setting could be made part of the global setting
    % <SETTINGS_OPTIM>. However, I feared that that would make the code
    % unreadable (these settings only make sense in their context), and is
    % not really necessary. However, for full flexibility, it can be
    % considered to put them in the global.
    
    if ind_final >= n_conf(1) % do we already have enough parameter sets in the total joint CI?
        if ind_single >= n_conf(2) % and also enough in the inner rim?
            flag_stop = 1; % then we can stop!
            
        else % what to do if there are enough within the outer rim, but few in the inner?
            disp('  Next round will focus on inner rim (outer rim has enough points)')
            coll_ok = coll_all(1:ind_inner,:); % only take the ones that are within the inner rim (plus a little extra)
            flag_inner = 1; % signal that we'll do inner rim only in next round
            if ind_inner < n_ok % if there are very few currently in the inner rim, take the <n_ok> best ones from <coll_all> ...
                coll_ok = coll_all(1:n_ok,:); % take the best <n_ok> values to continue to next round
            elseif ind_inner > 0.5 * n_conf(2) % if we already have quite some values in inner rim ...
                % Focus on the <n_ok> sets close to the rim (now, <n_ok> on both sides of the cut-off)
                coll_ok = coll_all(max(1,ind_single-n_ok):ind_single+n_ok,:);
                coll_ok  = cat(1,coll_ok,coll_all(1,:)); % also add the optimised best value in there ... (helps to look for better optimum)
            end
        end
        
    else % then there are not enough accepted parameter sets in the total joint CI
        
        if ind_cont_t > n_ok % if there are enough in <coll_tries> to continue with ...
            % in principle, the next round will continue from sets in
            % <coll_tries> and not from <coll_all>. The reason is that it is
            % better to propagate new sets than to continue propagating
            % sets that have been propagated before already.
            
            if ind_cont_t > 2 * n_conf(1) % if we will try more than 2 times of what we finally need ...
                
                if ind_final > 0.5 * n_conf(1) % if we already have half the values we need
                    % focus on values around the edge of the inner rim as
                    % that is where we need most precision.
                    coll_ok = coll_tries([max(1,ind_single-n_ok):min(ind_cont_t,ind_single+n_ok)],:);
                    
                else % then we have a lot of values to try, but not so much accepted yet
                    
                    if coll_tries(2*n_conf(1),end) > mll + chicrit_max % at 2 * <n_conf(1)>, do we still include some bad points?
                        ind_cont_t = 2 * n_conf(1); % limit continuation points to the best 2 * <n_conf(1)>
                    else % then it is not a good idea to limit continuation to 2 * <n_conf(1)>
                        ind_cont_t  = ind_cont_t2; % take the sets within <chicrit_max>
                        % This could still be a lot ... but then we'll
                        % trigger a reduction in new mutations.
                    end
                    coll_ok  = coll_tries(1:ind_cont_t,:);
                end
            else % then we don't try too much in the next round, so just use <ind_cont_t> as planned
                coll_ok  = coll_tries(1:ind_cont_t,:);
            end
            coll_ok  = cat(1,coll_ok,coll_all(1,:)); % also add the optimised best value in there ... (helps to look for better optimum)
            
        else % there are NOT enough values in <coll_tries> to continu with ... then we need to look at <coll_all>.
            if ind_cont_a > n_ok % if there are enough reasonable sets in total <coll_all>
                if ind_cont_a > 2 * n_conf(1) % if we will try more than two times of what we finally need ...
                    if coll_all(2 * n_conf(1),end) > mll + chicrit_max % at 2 * <n_conf(1)>, do we still include some bad points?
                        ind_cont_a = 2 * n_conf(1); % limit to best 2 * <n_conf(1)> to continue
                    else
                        ind_cont_a = ind_cont_a2; % otherwise take the ones within <crit_max>
                        % this could still be a lot ...
                    end
                end
                coll_ok  = coll_all(1:ind_cont_a,:); % continue with selection from the total <coll_all>
            else % if all else fails ...
                coll_ok  = coll_all(1:n_ok,:); % take the best <n_ok> values from the total <coll_all>
            end
        end
    end
    
    % Now that we know the set that we will continue with, it is time to
    % update the plot of parameter space. And make a plot of the progress
    % so far
    if plot_intermed == 1
        if flag_inner == 0 && flag_stop ~= 1
            figh = plot_grid(pmat,coll_ok,coll_all,[],figh,SETTINGS_OPTIM);
        else % if we were refining the inner rim, or will stop, we can plot the outer rim instead
            figh = plot_grid(pmat,[],coll_all,[],figh,SETTINGS_OPTIM);
        end
        uistack(f,'top')  % but place progress bar on top!
    end
    
    % BLOCK 4.6. See if we need a new round, and if so, prepare for it.
    % Again, special care is taken to avoid taking too few or too many new
    % tries in the next round. Therefore, the tricky part is to come up
    % with an efficient value for <n_tr_i> (number of random mutations per
    % set in the next round).
    %
    % Note: the settings here are all tweaked to get good results in the
    % cases tested. These setting could be made part of the global setting
    % <SETTINGS_OPTIM>. However, I feared that that would make the code
    % unreadable (these settings only make sense in their context), and is
    % not really necessary. However, for full flexibility, it can be
    % considered to put them in the global.
    
    if n_rnd == n_max && flag_stop == 0 % at some point, we need to force a stop ... and examine what went wrong
        flag_stop = 1; % let's stop here
        disp(' ')
        disp(['We have now done ',num2str(n_rnd),' rounds without reaching the stopping criterion ... we stop here!'])
    end
        
    % If we need another round, prepare for it.
    if flag_stop ~= 1
        n_rnd = n_rnd + 1; % increase counter for rounds by 1
        if n_rnd <= length(n_tr)     % do we still have values for the optim. criteria left in our options set?
            n_tr_i    = n_tr(n_rnd); % number of mutations per set
            f_d_i     = f_d(n_rnd);  % maximum step as factor of grid spacing for mutations
            chicrit_i = chicrit_rnd(n_rnd); % criterion to select ok values for continuation
        else % otherwise, take the last values for any additional round
            n_tr_i    = n_tr(end);
            f_d_i     = f_d(end);
            chicrit_i = chicrit_rnd(end);
        end
        
        % Modify number of tries if we have a lot or only few sets within the total cloud or inner rim.
        crit_ntry = [ind_final/n_conf(1) ind_single/n_conf(2)]; % vector: how far are we from target number of sets in total and inner rim?
        if crit_ntry(1+flag_inner) > 0.75 || size(coll_ok,1) > 2000
            % If we already have more than 75% of what we finally need, or more than 2000 sets to continue with
            % Note: when <flag_inner> = 1, we look at the status of the inner rim!
            n_tr_i = floor(n_tr_i/2); % decrease the number of tries per set for next round
        elseif n_rnd > 3 % starting at round 4, we start to worry if we haven't found so many yet
            if size(coll_ok,1) < 0.5 * n_conf(1+flag_inner) && size(coll_ok,1) < 1000
                % If we have less than 50% of what we finally need, and less than 1000 sets to try next ...
                n_tr_i = 2 * n_tr_i; % use twice as many tries in the next round
                if n_rnd > 4 && crit_ntry(1+flag_inner) < 0.25 % if there is very few when we start round 5 ...
                    n_tr_i = 2 * n_tr_i; % use twice as many tries (again) in the next round
                end
            elseif n_rnd > 7 % in the last few rounds, always worry if we still have less than 75% of the final number ...
                n_tr_i = 2 * n_tr_i; % use twice as many tries in the next round
            end
        end
        n_tr_i = min(n_tr_i,max(2,floor(10*n_conf(1+flag_inner)/size(coll_ok,1)))); % limit the number of tries per set for next round
        % This allows a max of 10x the target value to be tried in the next
        % round (scaling back <n_tr_i>). Well, sometimes more, as at least
        % 2 tries will be allowed for each set in <coll_ok>. I used to have
        % a check for whether we are in the final rounds or if we found
        % very little, but I don't think that is needed (and sometimes
        % leads to large numbers of new tries). If all fails, the extra
        % sampling rounds will come to the rescue.
        
    end
    
    % Prune <coll_all> to remove all values that are outside highest criterion
    coll_all = coll_all(coll_all(:,end) < mll + chicrit_max,:);
   
end

% BLOCK 4.7. We only get here when the while loop is finished. To improve
% upon the best fit ... do a thorough optimisation

disp('  Finished sampling, running a simplex optimisation ...')

pfit = coll_all(1,1:end-1)'; % copy best values to <pfit>, for parameters that are to be fitted
% pmat(ind_fit,1) = pfit; % update <pmat> with the best value at this point NOT NEEDED
[phat,mll]      = setup_simplex(pfit,1,pmat); % do a normal optimisation
pmat(ind_fit,1) = phat; % update <pmat> with the best fit parameters
coll_all        = [phat' mll; coll_all]; % place the new optimum at the top of <coll_all>

% Update the location of the total and inner conf. region.
ind_final  = find(coll_all(:,end) < mll + chicrit_joint,1,'last');
ind_single = find(coll_all(:,end) < mll + chicrit_single,1,'last');

disp(['  Status: ',num2str(ind_final),' sets within total CI and ',num2str(ind_single),' within inner. Best fit: ',num2str(coll_all(1,end))])

%% BLOCK 4B. Specific for BYOM!
% We can stop here, if only a rough search is needed ... then we skip
% profiling and additional sampling, and only use the basic sampling
% rounds. This is useful for slow analyses (e.g., DEBtox) and if we are
% mainly interested in the best fit and a rough idea of the parameter
% landscape. For proper analysis, stopping here is NOT recommended, since
% there are cases where the initial sampling misses a relevant part of
% parameter space.

if opt_optim.ps_profs == 0 % only if option says we can skip profiling
    
    % Calculate how many sets will be used for propagation.
    ind_prop1 = find(coll_all(:,end) < coll_all(1,end) + 0.5 * SETTINGS_OPTIM.crit_prop(1),1,'last');
    ind_prop2 = find(coll_all(:,end) < coll_all(1,end) + 0.5 * SETTINGS_OPTIM.crit_prop(2),1,'last');
    if isempty(ind_prop1) % <SETTINGS_OPTIM.crit_prop(1)> could be zero (if we want entire inner rim), which leads to empty <ind_prop1>
        ind_prop1 = 1; % start from the first value
    end
    
    % Start a matrix <pmat_print> that will collect the best parameters and CIs
    % (and additional info for each parameter).
    pmat_print = (coll_all(1,1:end-1) )'; % add best values so far in first column
    if ind_prop2 > ind_prop1
        pmat_print(:,[2 3]) = [(min(coll_all(ind_prop1:ind_prop2,1:end-1)))' (max(coll_all(ind_prop1:ind_prop2,1:end-1)))']; % add the edges of the propagation cloud (approx. single-parameter CIs)
    else % there may be cases where all elements of coll_all are already within ind_prop1
        pmat_print(:,[2 3]) = NaN; % just use NaN: CIs cannot be estimated
    end
    pmat_print(:,[4 5]) = [(min(coll_all(1:ind_final,1:end-1)))' (max(coll_all(1:ind_final,1:end-1)))']; % add the edges of the joint 95% cloud
    pmat_print(:,[6 7 8]) = zeros(n_fit,3); % add three columns with zeros (used to indicate where bounds are hit and broken CIs)
    pmat_print(ind_log,1:5) = 10.^(pmat_print(ind_log,1:5)); % put log parameters back on normal scale (also their CIs and total bounds)
    % NOTE: UNLIKE PMAT THIS WILL THUS BE ON NORMAL SCALE

    for i_p = 1:length(ind_fit) % trick to make sure profiles are plotted (but without refinement red line)
        coll_prof_pruned{i_p} = [NaN NaN];
    end
    
    % And make a plot of the final results as parameter-space plot
    figh = plot_grid(pmat,[],coll_all,coll_prof_pruned,figh,SETTINGS_OPTIM); % should be superfluous
    
    figure(figh)  % use existing handle and make sure that the multiplot is the current plot
    % Add a title to the plot.
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    text(0.5, 1,['File: ',glo.basenm,' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
    if glo.saveplt > 0
        save_plot(figh,['parspace_',glo.basenm]) % save parspace plot in output folder
    end
    snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output
    
    % Save the sample into a MAT file, so it can be used to create CIs on
    % model curves. Prune <coll_all> to remove all values that are outside
    % highest criterion.
    coll_all = coll_all(coll_all(:,end) < mll + max(chicrit_joint,0.5*SETTINGS_OPTIM.crit_prop(2)),:);
    
    % Note that we now also save a reconstructed par with the mat file. First,
    % put parameters that need to be fitted on log scale back on normal scale
    pmat_tmp = pmat; % don't mess with the one that will be saved and output of the function!
    pmat_tmp(pmat_tmp(:,5)==0,1) = 10.^(pmat_tmp(pmat_tmp(:,5)==0,1));
    par = packunpack(2,[],pmat_tmp); % pack pmat_tmp into structure par
    
    save([glo.basenm,'_PS'],'pmat','par','coll_all','pmat_print','coll_prof_pruned','names')
    
    % If we're doing a DEBtox analysis, it is a good idea to save the
    % brood-pouch delay with the sample. 
    if isfield(glo,'Tbp')
        Tbp = glo.Tbp;
        save([glo.basenm,'_PS'],'Tbp','-append')
    end
    % If we have multiple versions of the same parameter, also save that info.
    if isfield(glo,'names_sep')
        names_sep = glo.names_sep;
        save([glo.basenm,'_PS'],'names_sep','-append')
    end
    
    stats = [ind_final,ind_single,1+ind_prop2 - ind_prop1]; % useful statistics for <calc_optim> to print
    delete(f) % delete progress bar again
    
    % create a pmat version that can be copied-pasted into the code for
    % refining the parspace run later
    pmat_disp = pmat;
    pmat_disp(ind_fit,1) = (coll_all(1,1:end-1) )'; % add best values so far in first column
    pmat_disp(ind_log_all,1) = 10.^pmat_disp(ind_log_all,1); % put all log-scale parameters back to normal scale
    % refine the search ranges (add 20% to the 95% joint cloud, but stay within original range)
    pmat_disp(ind_fit,3) = max(pmat_print(:,4)/1.2,pmat(ind_fit,3));
    pmat_disp(ind_fit,4) = min(pmat_print(:,5)*1.2,pmat(ind_fit,4));
    
    disp(' ')
    disp('Results from the limited parspace exploration (initial sampling only).')
    disp('This can be copied-pasted into your BYOM script.')
    disp('Warning: for oddly-shaped parameter spaces, this may miss local minima!')
    disp('(in extreme cases, even the global optimum can be missed)')
    disp('Suggestions for min-max search bounds returned in par_out.')
    disp('=======================================================================')
    print_par(pmat_disp,-2) % print in formatted manner, but only fitted parameters
    disp('=======================================================================')
    
    pmat(:,[3 4]) = pmat_disp(:,[3 4]); % this will return the new bounds in pmat, so we can automatically re-run!
    
    return % and go back to where we're called from
    
end


%% BLOCK 5. Calculation of profile likelihood.
% From the sample, calculate a profile likelihood as a line in a plot of
% MLL versus a single parameter. The line is the lowest value of the MLL
% that can be reached at each value of the parameter (target parameter
% fixed, remaining parameters fitted).

n_rnd = n_rnd + 1; % increase counter for rounds by 1
disp(['Starting round ',num2str(n_rnd),', creating the profile likelihoods for each parameter'])
% Note: waitbar updating is dealt with in <calc_proflik_ps> now

% Call function to do the profiling. If profiling finds a lower MLL, <phat>
% returns the new best-fitting parameters, and <mll> the new MLL (otherwise,
% the previous MLL is returned). <coll_prof> will be the profile info for
% each parameter (cell array), and <coll_ok> the parameter sets that are
% candidates for further sampling.
[coll_prof,coll_ok,phat,mll,coll_prof_pruned,flag_short] = calc_proflik_ps(pmat,coll_all,bnds_tmp,f,SETTINGS_OPTIM);

if ~isempty(phat) % if profiling led to a better optimum
    pmat(ind_fit,1)= phat; % update <pmat> with the better fitting parameters
    coll_all = [pmat(ind_fit,1)' mll; coll_all]; % place the new optimum at the top of <coll_all>
    if coll_all(2,end) - coll_all(1,end) > SETTINGS_OPTIM.real_better % if there is a real improvement
        coll_ok  = [coll_all(1,:) ; coll_ok]; % also add the new optimum to <coll_ok>, so trigger resampling if <coll_ok> was empty
    end
    % Update the location of the total and inner conf. region.
    ind_final  = find(coll_all(:,end) < mll + chicrit_joint,1,'last');
    ind_single = find(coll_all(:,end) < mll + chicrit_single,1,'last');
end

% Code below is probably not needed anymore: ind_final may be 1, but
% <coll_ok> may still include useful sets. See also BLOCK 6.1.
% % NOTE: things go wrong if ind_final = 1, which may occur in very extreme
% % situations. For now, just produce an error and print new optimum on screen.
% if ind_final == 1
%     phat = coll_all(1,1:end-1);
%     pmat(ind_fit,1) = phat';
%     pmat(ind_log,1) = 10.^(pmat(ind_log,1));
%     disp('Parameters for the new optimum')
%     disp('=====================================')
%     for i = 1:length(ind_fit)
%         disp([glo2.names{ind_fit(i)},' ',num2str(pmat(ind_fit(i),1))])
%     end
%     disp('=====================================')
%     error('Cannot continue as there is only one element of the sample left to work with. Restart with tighter ranges.')
% end

if plot_intermed == 1
    figh = plot_grid(pmat,[],coll_all,coll_prof_pruned,figh,SETTINGS_OPTIM); % update the plot of parameter space, now with profiles
    % % TEMP save preliminary plot, just after first profiling (for design doc)
    % figure(figh) % use existing handle and make sure that the multiplot is the current plot
    % save_plot(figh,['output_report',filesep,[savestr,'_pre']],3) % save parspace plot as pdf in output folder
    uistack(f,'top')  % but place progress bar on top!
end
disp(['  Status: ',num2str(ind_final),' sets within total CI and ',num2str(ind_single),' within inner. Best fit: ',num2str(coll_all(1,end))])

%% BLOCK 6. Do an extra round of sampling when profiling indicates problems.
% The profiling will identify points where there is a difference between
% the profile and the sampling points. This generally indicates that
% sampling was not sufficiently covering parameter space. These points are
% returned in <coll_ok>, which we can use for additional rounds of sampling.

if ~isempty(coll_ok) || ind_single < n_conf(2) % if we have candidate sets for additional sampling, or if we don't have enough in inner rim
    % When profiling finds a much better optimum, we can suddenly end up
    % with insufficient sets in the inner rim (as it will shift the profile
    % up). 
    
    % BLOCK 6.1. Initial things when we start the extra sampling.
    
    if isempty(coll_ok) || (ind_single < n_conf(2) && size(coll_ok,1) < 10) % when empty or when there are too few in inner rim (happens when we found a much better optimum in profiling)
        % Add values around the edge of the (new) inner rim as that is where we need most precision.
        coll_ok = [coll_ok ; coll_all(max(1,ind_single-n_ok):min(ind_single+n_ok,size(coll_all,1)),:)];
        % This is now modified with a <min> statement to catch cases where
        % <coll_all> does not have the required number of elements.
    end
    
    % NOTE: in some extreme cases, coll_all may be very small, so also
    % include edges of coll_ok to find the edges! (that includes the points
    % flagged by profiling as well). This will not help when coll_ok is
    % empty, so there may still be an error at some point ...
    coll_tmp    = [coll_ok(:,1:end-1) ; coll_all(1:ind_final,1:end-1)]; % this is the matrix to derive new cloud edges from
    edges_cloud = [(min(coll_tmp,[],1))' (min(coll_tmp,[],1))']; % edges of the joint 95% cloud
    clear coll_tmp % clear the temporary variable for memory use
    
    d_grid      = SETTINGS_OPTIM.d_extra * diff(edges_cloud,1,2); % create a new grid spacing as fraction of the new total parameter range
    f_d_i       = SETTINGS_OPTIM.f_d_extra; % maximum jump size as fraction of the grid spacing
    
    coll_ok     = sortrows(coll_ok,n_fit+1); % sort <coll_ok> based on minloglik
    coll_ok     = coll_ok(coll_ok(:,end) < mll + SETTINGS_OPTIM.crit_add_extra,:); % prune <coll_ok> to remove all values that are too far outside inner rim
    n_rnd_x     = 0; % set simple counter for extra rounds to zero
    n_test      = 0; % a simple counter for how many times we trigger new profiling
    
    if ind_single < n_conf(2)
        disp('  Profiling has led to a new optimum, which left insufficient sets in the inner rim: extra sampling rounds will be started.')
    elseif ~isempty(coll_ok) % in some exceptional cases, <coll_ok> is left empty after the pruning step above
        disp('  Profiling has detected gaps between profile and sample, which requires extra sampling rounds.')
    end
    
    % BLOCK 6.2. Continue sampling and testing until there is nothing left
    % to worry about (and <coll_ok> will be empty).
    
    while ~isempty(coll_ok) % keep going until we are not given new sets to mutate
        
        % BLOCK 6.2.1. Prepare for a new round of sampling.
        n_cont  = size(coll_ok,1); % how many sets to continue with
        n_rnd   = n_rnd + 1;   % increase counter for rounds by 1
        n_rnd_x = n_rnd_x + 1; % increase counter for extra rounds by 1
        
        % Select number of mutations per parameter set in <coll_ok>.
        n_tr_i  = 40; % basic number of new tries per set
        n_tr_i  = min(n_tr_i,max(2,floor(5*n_conf(1) / n_cont))); % make sure that we don't have more than a fixed number new tries
        % This scales <n_tr_i> back if needed (avoids extreme amounts of
        % mutations) to a max. of 5 times the value we aim for in the total
        % cloud (with a minimum of 2 tries per set).
        
        disp(['Starting round ',num2str(n_rnd),', (extra ',num2str(n_rnd_x),') refining a selection of ',num2str(size(coll_ok,1)),' parameter sets, with ',num2str(n_tr_i),' tries each'])

        % BLOCK 6.2.2. Call <rand_mutations> to mutate <coll_ok>, give each
        % new set an MLL, and return it in <coll_tries>.
        coll_tries = rand_mutations(pmat,coll_ok,bnds_tmp,n_tr_i,f_d_i*d_grid,f,n_rnd);

        coll_all   = cat(1,coll_all,coll_tries); % add the tries to <coll_all>
        coll_all   = sortrows(coll_all,n_fit+1); % also sort based on minloglik
        
        % Prune <coll_all> to remove all values that are outside highest criterion
        coll_all = coll_all(coll_all(:,end) < mll + chicrit_max,:);
        
        if plot_intermed == 1
            figh = plot_grid(pmat,[],coll_all,coll_prof_pruned,figh,SETTINGS_OPTIM); % update plot of parameter space
            uistack(f,'top')  % but place progress bar on top!
        end
        
        % Update the location of the total and inner conf. region.
        ind_final  = find(coll_all(:,end) < mll + chicrit_joint,1,'last');
        ind_single = find(coll_all(:,end) < mll + chicrit_single,1,'last');
        
        disp(['  Status: ',num2str(ind_final),' sets within total CI and ',num2str(ind_single),' within inner. Best fit: ',num2str(coll_all(1,end))])
        
        % BLOCK 6.2.3. Run a test to see whether the sample now DOES match
        % the profile likelihood. It checks where the profile is
        % considerably lower than the sample (which means that more
        % sampling is needed), and when the sample is lower than the
        % profile (which means we may need a new profile).
        [flag_profile,coll_ok] = test_proflik_ps(pmat,coll_prof,coll_all,SETTINGS_OPTIM);
         
        % BLOCK 6.2.4. Do we need another round of extra sampling? First,
        % make sure that there are enough points in the inner rim
        % (otherwise, continue sampling).
        if ind_single < n_conf(2) % if we don't have enough in inner rim ... simply continue sampling until we do
            % Theoretically, this could lead to an infinite loop. However,
            % I don't think that there can be cases where continuous
            % sampling won't fill the inner rim.
            ind_cont = find(coll_tries(:,end) < mll + SETTINGS_OPTIM.crit_add_extra,1,'last'); % only continue with new tries that are close to inner rim
            if ind_cont < 0.5 * n_conf(2) % if tries are less than half of what we finally need ...
                coll_ok = [coll_ok ; coll_tries(1:ind_cont,:)]; % also add the <coll_ok> to it
            else
                coll_ok = coll_tries(1:ind_cont,:); % otherwise, only tries
            end
            f_d_i = max(0.1,f_d_i * 0.7); % decrease jumps by a factor over subsequent rounds to let extra sampling contract
            
        elseif ~isempty(coll_ok) && n_rnd_x < 5
            % BLOCK 6.2.5. We already have enough points in the inner rim,
            % but gaps have been flagged and we haven't done more than 4
            % extra sampling rounds yet. If there are still gaps left after
            % 5 rounds, it is unlikely that more sampling will really help.
            
            % Continue sampling with <coll_ok>.
            f_d_i = max(0.1,f_d_i * 0.7); % decrease jumps by a factor over subsequent rounds to let extra sampling contract
        
        elseif (flag_profile(2) > 0.05 || flag_short == 1) && n_test < 2 % if some difference between profile and sample is spotted
            % BLOCK 6.2.6. At this point, see if we have a flag for a new
            % profile. The profile flag is a vector of 2: number of points
            % flagged, and maximum distance between profile and sample. At
            % this point, I only trigger on the latter (note that I allow
            % some deviation). The inclusion of <n_test> makes sure that we
            % don't end up in an infinite loop ... if we already did 2
            % reprofiles, another one will probably not help. Also do a new
            % profiling if the previous one indicated that one or more
            % profiles were short.
            %
            % Note: there are some data sets that cannot be profiled
            % properly with the current algorithm. They continue to trigger
            % reprofiling. The consequences for the user are minor: a small
            % error on the single parameter CI and a bit longer calibration
            % time (and a sample that is a bit larger than needed).
            % Therefore, I think there is no pressing need to add more code
            % to catch these cases.
            
            n_rnd   = n_rnd + 1; % increase counter for rounds by 1
            disp(['Starting round ',num2str(n_rnd),', refining the profile likelihoods for each parameter again'])
            n_test = n_test + 1; % increase the counter for number of re-profiles
            
            % Call function to do the profiling. If profiling finds a
            % lower MLL, <phat> returns the new best-fitting parameters
            % and <mll> the new MLL (otherwise, the previous MLL is
            % returned).
            [coll_prof,coll_ok,phat,mll,coll_prof_pruned,flag_short] = calc_proflik_ps(pmat,coll_all,bnds_tmp,f,SETTINGS_OPTIM);

            if ~isempty(phat) % if profiling led to a better optimum
                pmat(ind_fit,1)= phat; % update <pmat> with the better fitting parameters
                coll_all = [pmat(ind_fit,1)' mll; coll_all]; % place the new optimum at the top of <coll_all>
                % Update the location of the total and inner conf. region.
                ind_final  = find(coll_all(:,end) < mll + chicrit_joint,1,'last');
                ind_single = find(coll_all(:,end) < mll + chicrit_single,1,'last');
            end
            
            disp(['  Status: ',num2str(ind_final),' sets within total CI and ',num2str(ind_single),' within inner. Best fit: ',num2str(coll_all(1,end))])
            if plot_intermed == 1
                figh = plot_grid(pmat,[],coll_all,coll_prof_pruned,figh,SETTINGS_OPTIM); % update plot of parameter space with profile
                uistack(f,'top')  % but place progress bar on top!
            end
            
            % If there is something in <coll_ok>, found by profiling, a
            % new round will be triggered. Then we need to reset the
            % extra rounds counter and jump size.
            if ~isempty(coll_ok)
                n_rnd_x = 0; % reset extra rounds counter
                f_d_i   = SETTINGS_OPTIM.f_d_extra; % reset maximum jump size as fraction of the grid spacing
                disp('  New profiling has detected gaps, which requires additional extra sampling rounds!')
            end
            
        else
            % BLOCK 6.2.7. If nothing is flagged, we can stop, but display
            % on screen if we stopped at maximum number of rounds.
            if size(coll_ok,1) > 0 
                disp(['Exiting parameter space explorer, but still ',num2str(size(coll_ok,1)),' sets flagged (too much distance between profile and sample). Check plot.'])
            end
            if flag_profile(2) > 0.05
                disp('Exiting parameter space explorer, but there still seem to be parameter sets outside of the profile curves. Check plot.')
            end
            coll_ok = []; % make <coll_ok> empty (if it is not empty already) so we can stop further sampling
        end
        
    end
end

% Prune <coll_all> to remove all values that are outside highest criterion.
% This is needed when a (much) better optimum has been found by extra
% sampling.
coll_all = coll_all(coll_all(:,end) < mll + chicrit_max,:);

%% BLOCK 7. Calculate CIs for single parameters
% Now we can extract the CIs for the single parameters. Do this separately
% from the profiling because the profiling (and additional sampling) might
% just find a better optimum (which shifts all profiles, and thereby the
% CIs).

% Start a matrix <pmat_print> that will collect the best parameters and CIs
% (and additional info for each parameter).
pmat_print = (coll_all(1,1:end-1) )'; % add best values so far in first column
pmat_print(:,[2 3]) = nan; % add two columns for the single-parameter CIs
pmat_print(:,[4 5]) = [(min(coll_all(1:ind_final,1:end-1)))' (max(coll_all(1:ind_final,1:end-1)))']; % add the edges of the joint 95% cloud
pmat_print(:,[6 7 8]) = zeros(n_fit,3); % add three columns with zeros (used to indicate where bounds are hit and broken CIs)

for i_p = 1:n_fit % run through all fitted parameters
    
    parprof  = coll_prof{i_p,1}; % read <parprof> again from the collected ones
    coll_all = cat(1,coll_all,parprof); % add the profile information to the sample in <coll_all>!
    coll_all = sortrows(coll_all,n_fit+1); % sort <coll_all> again based on minloglik, as points have been added
    % This may help if the sampling was not very good. I do this here now,
    % as otherwise we cannot do the profiling a second time (the first
    % profile will be included in the sample, so that is always the best
    % sample point in each slice).
    
    % BLOCK 7.1. Check if *other* parameters are running into bounds at the
    % min/max of *this* parameter (that could have influenced the bounds).
    % This is not perfect, though it will catch the obvious cases.
    %
    % THIS BLOCK HAS BEEN REMOVED: IT DID NOT WORK WELL ENOUGH, SO WE'LL
    % LEAVE THAT FOR FUTURE VERSIONS IF NEEDED.
    
    % BLOCK 7.2. Find the outer bounds of the CI (we don't care for broken
    % intervals anymore; they are flagged, not calculated). The CI bounds
    % themselves are collected in <pmat_print> (column 2 and 3). Column 6
    % and 7 will get a 1 when the CI is open on the lower or upper end.
    % Column 8 will get a 1 when there are multiple zero-crossings and we
    % have a broken CI.
    %
    % Note: it would also be possible to use BYOM's <calc_xing> for this
    % purpose (and report all crossings). However, for now, it is probably
    % better to stick close to the openGUTS version.
    
    prof_tst = coll_prof_pruned{i_p}; % use the pruned version to calculate CIs
    prof_tst(:,2) = prof_tst(:,2)-mll-chicrit_single; % rephrase the problem to zero finding
    
    prof_tst(~isfinite(prof_tst(:,2)),:) = []; % remove any rows where the MLL is inf or nan
    % Note: INF should not happen anymore (pruned in <prune_prof>), but
    % perhaps NaNs could happen?
    
    % First and last value in <prof_tst> should be positive, otherwise we
    % have an open interval.
    if prof_tst(1,2) < 0 % that means that the CI is open on the lower end: parameter is hitting lower bound
        pmat_print(i_p,2) = bnds_tmp(i_p,1); % just use lower bound of parameter range
        pmat_print(i_p,6) = 1; % and mark it in <pmat_print>
    else % otherwise, we need to find the first crossing point
        ind_low = find(prof_tst(:,2)<0,1,'first'); % find index to first point below cut_off
        % Next, linearly interpolate between points just above and below cut-off.
        pmat_print(i_p,2) = interp1([prof_tst(ind_low-1,2) prof_tst(ind_low,2)],[prof_tst(ind_low-1,1) prof_tst(ind_low,1)],0);
    end
    if prof_tst(end,2) < 0 % that means that the CI is open on the upper end: parameter is hitting upper bound
        pmat_print(i_p,3) = bnds_tmp(i_p,2); % just use upper bound of parameter range
        pmat_print(i_p,7) = 1; % and mark it in <pmat_print>
    else
        ind_hi = find(prof_tst(:,2)<0,1,'last'); % find index to last point below cut_off
        % Next, interpolate between points just above and below cut-off.
        pmat_print(i_p,3) = interp1([prof_tst(ind_hi,2) prof_tst(ind_hi+1,2)],[prof_tst(ind_hi,1) prof_tst(ind_hi+1,1)],0);
    end
    
    % And a final count of the number of zero crossings. When there are two
    % we have nice regular CI.
    num_x = sum(diff(sign(prof_tst(:,2)))~=0); % this counts the number of sign changes in the profile to test
    num_x = num_x + sum(pmat_print(i_p,[6 7]) == 1); % also add number of open intervals as a zero crossing
    if num_x > 2 % then we have a broken CI
        pmat_print(i_p,8) = 1; % mark it in <pmat_print>
    end
    
end

%% BLOCK 8. Final plotting saving the results to a MAT file.
% save the sample in a MAT-file with the name of the script, with _SD or
% _IT at the end. <load_sample> can be used to retrieve it.

% Now we can safely put <pmat_print> on normal scale, where needed.
pmat_print(ind_log,1:5) = 10.^(pmat_print(ind_log,1:5)); % put log parameters back on normal scale (also their CIs and total bounds)
% NOTE: UNLIKE PMAT THIS WILL THUS BE ON NORMAL SCALE

% And make a plot of the final results.
figh = plot_grid(pmat,[],coll_all,coll_prof_pruned,figh,SETTINGS_OPTIM); % should be superfluous
drawnow % make sure that everything is plotted that needs to be plotted

if opt_optim.ps_notitle == 0 % add a title to the plot.
    figure(figh) % use existing handle and make sure that the multiplot is the current plot    
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    text(0.5, 1,['File: ',glo.basenm,' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
end

if glo.saveplt > 0
    save_plot(figh,['parspace_',glo.basenm]) % save parspace plot in output folder
end
snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output

% Save the sample into a MAT file, so it can be used to create CIs on model
% curves.
% 
% Note that the first column of <pmat> stays on log-scale where needed.
% This is different from the procedure in openGUTS.

% Prune <coll_all> to remove all values that are outside highest criterion,
% and only keep the columns that are fitted (<load_sample> will then be
% able to reconstruct the entire matrix). We could decide to save only a
% partial sample here. For example, only the values that will be propagated
% (that makes the job of <load_sample> a bit easier. 
coll_all = coll_all(coll_all(:,end) < mll + max(chicrit_joint,0.5*SETTINGS_OPTIM.crit_prop(2)),:);
% The max function here will only come into play when we have one fitted
% parameter (otherwise, <chicrit_joint> is highest). Without it, for one
% fitted parameter, the propagation sample would be everything in the inner
% rim (which equals the joint CI for df=1), and not the part within the
% horizontal dotted lines.

% Note that we will save the pruned version of the profile likelihood
% <coll_prof_pruned>. At this moment, <coll_prof_pruned> is a cell array:
% each cell is a two-column matrix with parameter values and MLLs. This can
% be turned into a regular matrix if we like to save this information in a
% text file (but note that <coll_prof_pruned> can be of different length
% for different parameters).

names = glo2.names; % also save the <names> variable with the sample
ind_fittag   = ~strcmp(names,'tag_fitted');
names        = names(ind_fittag); % make sure that the fit tag is not in names_saved

% Note that we now also save a reconstructed par with the mat file. First,
% put parameters that need to be fitted on log scale back on normal scale
pmat_tmp = pmat; % don't mess with the one that will be saved and output of the function!
pmat_tmp(pmat_tmp(:,5)==0,1) = 10.^(pmat_tmp(pmat_tmp(:,5)==0,1));
par = packunpack(2,[],pmat_tmp); % pack pmat_tmp into structure par

save([glo.basenm,'_PS'],'pmat','par','coll_all','pmat_print','coll_prof_pruned','names')

% If we're doing a DEBtox analysis, it is a good idea to save the
% brood-pouch delay with the sample. 
if isfield(glo,'Tbp')
    Tbp = glo.Tbp;
    save([glo.basenm,'_PS'],'Tbp','-append')
end
% If we have multiple versions of the same parameter, also save that info.
if isfield(glo,'names_sep')
    names_sep = glo.names_sep;
    save([glo.basenm,'_PS'],'names_sep','-append')
end

% Also calculate how many sets will be used for propagation.
ind_prop1 = find(coll_all(:,end) < coll_all(1,end) + 0.5 * SETTINGS_OPTIM.crit_prop(1),1,'last'); 
ind_prop2 = find(coll_all(:,end) < coll_all(1,end) + 0.5 * SETTINGS_OPTIM.crit_prop(2),1,'last'); 
if isempty(ind_prop1) % <SETTINGS_OPTIM.crit_prop(1)> could be zero (if we want entire inner rim), which leads to empty <ind_prop1>
    ind_prop1 = 1; % start from the first value
end

% NEW after v.1.0: Update the location of the total and inner conf. region.
% Otherwise, <calc_optim> will report the status as achieved BEFORE adding
% the profile points to <coll_all>.
ind_final  = find(coll_all(:,end) < mll + chicrit_joint,1,'last');
ind_single = find(coll_all(:,end) < mll + chicrit_single,1,'last');

stats = [ind_final,ind_single,1+ind_prop2 - ind_prop1]; % useful statistics for <calc_optim> to print
delete(f) % delete progress bar again

