function [coll_prof,coll_ok,phat,mll,coll_prof_pruned,flag_short] = calc_proflik_ps(pmat,coll_all,bnds_tmp,f,SETTINGS_OPTIM)

% Usage: [coll_prof,coll_ok,phat,mll,coll_prof_pruned,flag_short] = calc_proflik_ps(pmat,coll_all,bnds_tmp,f,SETTINGS_OPTIM)
% 
% Use the collected sample from parameter space to zoom in on the profile
% likelihood for the single parameters. This is useful to extract
% single-parameter CIs, but also to extend the sample with profiled sets.
% The latter may be useful when sampling is poor in some part of parameter
% space. Furthermore, the profiling does a number of optimisations (one
% parameter fixed), which may spot a better value of the optimum.
%
% Strategy is to split up the parameter range for each parameter into 50
% slices of equal width. Within each slice, look for the best-fitting set.
% Use that set as starting point for an optimisation (with the target
% parameter fixed to the middle of the slice). Also do another optimisation
% from the optimised value in the previous slice. Sometimes, extending a
% profile works better than taking the best available sample in the slice.
% From these two optimisations, the best value is optimised again and
% stored. Next, the profiling starts from the other side and calulates the
% next version from the previous profile point. In some cases, extending a
% profile from the other direction can find a better profile.
%
% Additionally, the algorithm is looking for situations where there is a
% substantial difference between the profile and the sample. If that
% happens in a relevant range (where the fit is still acceptable), the
% points are collected and will be used for an additional round of sampling
% (apparently, sampling was not optimal in that case). Note that the entire
% profile will also be added to the total sample as they are relevant
% parameter sets (this happens in <calc_parspace> at the end, when making
% CIs).
% 
% Extra output is a pruned version of <coll_prof> that is more useful for
% plotting and saving.
% 
% Inputs
% <pmat>      parameter matrix
% <coll_all>  the entire sample (matrix)
% <bnds_tmp>  boundaries of parameter space (log-transformed were needed)
% <f>         handle for waiting bar
% <SETTINGS_OPTIM> = settings read from setup_settings by <calc_optim_ps>
% 
% Outputs
% <coll_prof>  profile likelihoods with all parameters for each point (cell array)
% <coll_ok>    matrix with all candidate parameter sets for resampling
% <phat>       if the profiling found a better overall optimum, this
%              collects it (vector with fitted parameters only)
% <mll>        this function also returns the MLL; usually, this is the top
%              value in <coll_all> (only different when a better optimum is
%              found)
% <coll_prof_pruned>  pruned profile information; cell array with 2-column
%                     matrix for each fitted parameter
% <flag_short>        flag that is 1 when relevant part of the profile of
%                     at least one parameter is too short (triggers
%                     reprofiling)
% 
% Author     : Tjalling Jager
% Date       : May 2020
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <calc_proflik> code that is
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

%% BLOCK 1. Initial things.

global glo2

names = glo2.names; % names of the model parameters

ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
n_fit   = sum(pmat(:,2));       % number of fitted parameters
gridsz  = 50;                   % number of points for the profile over the parameter's range

% Note: options for rough/detailed simplex optimisations are set in <setup_simplex>

mll     = coll_all(1,end); % read the best min log-likelihood from first position of <coll_all>
mll_rem = mll; % remember the initial one entered (<mll> may be updated during profiling)
chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1);     % criterion for single-parameter CIs
chicrit_joint  = 0.5 * SETTINGS_OPTIM.crit_table(n_fit,1); % criterion for joint 95% CI or parameters
chicrit_joint  = max(chicrit_joint,chicrit_single+1.5); % makes sure to take a bit more than the single-par cut-off for 1 and 2 fitted parameters

ind_final   = find(coll_all(:,end) < mll + chicrit_joint,1,'last'); % find index to last point still in 95%CI
coll_prof   = cell(n_fit,1); % initialise an empty cell array that will collect the profiles for plotting later
coll_ok     = []; % start with an empty matrix to collect candidate sets that will go through another round of sampling
par_best    = []; % collect best better parameter set when profiling finds a better optimum
edges_cloud = [(min(coll_all(1:ind_final,1:end-1)))' (max(coll_all(1:ind_final,1:end-1)))']; % edges of the joint 95% cloud (matrix)

%% BLOCK 2. Do the actual profiling.

for i_p = 1:n_fit % run through all fitted parameters
    
    % BLOCK 2.1. Initial things.
    parnr   = ind_fit(i_p); % index of this fitted parameter in <pmat>
    par_rng = linspace(edges_cloud(i_p,1),edges_cloud(i_p,2),gridsz); % vector of evenly-spaced values across the parameters total range (final cloud)
    gridsp  = 0.5 * diff(edges_cloud(i_p,:))/(gridsz-1); % half distance between points in the grid for this parameter
    
    % The profile needs to include the best value exactly, and we do this
    % by slightly shifting the range to the right.
    ind_best = find(par_rng<=coll_all(1,i_p),1,'last'); % find the value in the grid range just below the best value
    par_rng  = par_rng + (coll_all(1,i_p) - par_rng(ind_best)); % shift the range to include the best value
    
    % The variables below with a '1' at the end relate to the optimisations
    % from the best point in the sample. Variable with a '2' relate to
    % optimisations from the previous point on the profile.
    pmat_tst1          = pmat; % copy <pmat> to temporary matrix (which will be modified)
    pmat_tst1(parnr,2) = 0; % don't fit this parameter as we are profiling it
    pmat_tst2          = pmat_tst1; % create another copy
    ind_fit_tst        = find(pmat_tst1(:,2)==1); % new index to the parameters to be fitted (vector)
    
    % BLOCK 2.2. Run through the grid points and optimise.
    
    parprof = nan(gridsz,n_fit+1); % initialise an empty matrix to catch the profile (all parameters plus their likelihood)
        
    % The first thing to do is to run over the grid just made (based on the
    % sample) and optimise in each interval to find the best point (which
    % is the profile likelihood).
    for i_g = 1:gridsz % run through all grid points (the x-axis of the profile likelihood)
        
        waitbar(((i_p-1)*gridsz+i_g)/(n_fit*gridsz),f,['Profiling fitted parameters']) % update progress bar
        
        % BLOCK 2.2.1. Find best value from <coll_all> in this sub-range.
        % This parameter set will be used as starting values for
        % optimisation from the sample. Note that <coll_all> is sorted, so
        % the first within this range is automatically the best.
        ind_tst = find(coll_all(:,i_p)>(par_rng(i_g)-gridsp) & coll_all(:,i_p)<(par_rng(i_g)+gridsp),1,'first');
        
        if isempty(ind_tst) % it can be empty if there are no sample points in this sub-range
            [~,ind_tst] = min(abs(coll_all(:,i_p) - par_rng(i_g))); % just take the closest point available
            mll_compare = +inf; % and set MLL here to INF (so it will be flagged for a gap)
        else
            mll_compare = coll_all(ind_tst,end); % this is the best MLL from the sample in this slice
        end
        
        % BLOCK 2.2.2. Generate a profile point from this sample point.
        pmat_tst1(ind_fit,1) = coll_all(ind_tst,1:end-1)'; % put parameter values for this point from the sample in the temporary <pmat_tst1>
        pmat_tst1(parnr,1)   = par_rng(i_g); % force the value on the grid for the parameter to be profiled
        p_tst1               = pmat_tst1(ind_fit_tst,1); % parameter estimates to be fitted for this point
        [phat_tst1,mll_tst1] = setup_simplex(p_tst1,0,pmat_tst1); % do a rough optimisation
        
        % BLOCK 2.2.3. Generate a profile point from the previous profile point.
        if i_g > 1 % then we have a previous value on the profile line
            % Then also use the previous value to see if propagating that one works better!
            pmat_tst2(ind_fit,1) = parprof(i_g-1,1:end-1)'; % put parameter values for this point from the profile in the temporary <pmat_tst2>
            pmat_tst2(parnr,1)   = par_rng(i_g); % force the value on the current grid point for the parameter to be profiled
            p_tst2               = pmat_tst2(ind_fit_tst,1); % parameter estimates to be fitted for this point
            [phat_tst2,mll_tst2] = setup_simplex(p_tst2,0,pmat_tst2); % do a rough optimisation
        else
            mll_tst2 = inf; % otherwise, min-log-likelihood value from previously profiled point is INF
        end
        
        % BLOCK 2.2.4. Compare the two estimates for the new profile point.
        if mll_tst2 < mll_tst1 % if previously profiled leads to the better likelihood ...
            phat_tst = phat_tst2; % use run from the previous profiled point
            pmat_tst = pmat_tst2;
        else % otherwise ...
            phat_tst = phat_tst1; % use run from sample point
            pmat_tst = pmat_tst1;
        end
        
        % BLOCK 2.2.5. Do a normal optimisation on the best value to make
        % sure we have the optimum.
        [phat_tst,mll_tst]      = setup_simplex(phat_tst,1,pmat_tst); % do a normal optimisation
        pmat_tst(ind_fit_tst,1) = phat_tst; % place the best-fitted parameters in <pmat_tst>
        parprof(i_g,:)          = [pmat_tst(ind_fit,1)' mll_tst]; % and add the result to <parprof>
        
        % BLOCK 2.2.6. Check for gaps and better optima.
        if mll_tst < mll+chicrit_single+1 % only if we are not clearly above the parameter's CI ...
            % We don't care when there is a gap between profile and sample
            % points in parts of parameter space that results in bad fits
            % anyway.
            
            % Do some tests to see if there are candidates for additional sampling
            if mll_compare - mll_tst > SETTINGS_OPTIM.gap_extra
                % There is a point, but the optimised value is better:
                % collect both the parameter values for the point and
                % the profile into <coll_ok>! Note: <mll_compare> will be INF
                % when there is no set in this slice at all. That is fine:
                % closest sample point and profile point will be collected.
                coll_ok = [coll_ok ; coll_all(ind_tst,:) ; pmat_tst(ind_fit,1)' mll_tst];
            end
            if mll_tst < mll % then we found a better optimum, so collect it
                par_best = parprof(i_g,:);
                if mll_rem - mll_tst > SETTINGS_OPTIM.real_better && mll - mll_tst > SETTINGS_OPTIM.real_better
                    % Only print on screen when the new optimum is really
                    % better than the original one, and better than the one
                    % that is now the best.
                    disp(['  Better optimum found when profiling ',names{parnr}, ': ',num2str(mll_tst),' (best was ',num2str(mll),')'])
                end
                mll = mll_tst;
            end
        end
        
    end
    
    % BLOCK 2.3. Run through the profile in reverse to see if it can be improved.
    for i_g = gridsz-1:-1:1 % run through all grid points (the x-axis of the profile likelihood) in reverse!
        mll_tst1             = parprof(i_g,end); % MLL at this point of the profile
        pmat_tst2(ind_fit,1) = parprof(i_g+1,1:end-1)'; % put parameter values for the previous point from the profile in <pmat_tst2> (we run in reverse)
        pmat_tst2(parnr,1)   = parprof(i_g,i_p); % force the value on the current value for the grid for the parameter to be profiled
        
        p_tst2               = pmat_tst2(ind_fit_tst,1); % parameter estimates to be fitted for this point
        [phat_tst2,mll_tst2] = setup_simplex(p_tst2,1,pmat_tst2); % do a standard optimisation
        pmat_tst2(ind_fit_tst,1) = phat_tst2; % place the best-fitted parameters in <pmat_tst2>
        
        if mll_tst2 < mll_tst1 % profiling from right to left provided a better value
            parprof(i_g,:) = [pmat_tst2(ind_fit,1)' mll_tst2]; % replace this entry in <parprof>
            if mll_tst1 - mll_tst2 > SETTINGS_OPTIM.gap_extra
                % If the optimised value is really better: collect the profiled set into <coll_ok>!
                coll_ok = [coll_ok ; parprof(i_g,:)];
            end
            if mll_tst2 < mll % then we even found a better optimum, so collect it
                par_best = parprof(i_g,:);
                if mll_rem - mll_tst2 > SETTINGS_OPTIM.real_better && mll - mll_tst2 > SETTINGS_OPTIM.real_better
                    % Only print on screen when the new optimum is really
                    % better than the original one, and better than the one
                    % that is now the best.
                    disp(['  Better optimum found when reverse profiling ',names{parnr}, ': ',num2str(mll_tst2),' (best was ',num2str(mll),')'])
                end
                mll = mll_tst2;
            end
        end
    end
    
    % BLOCK 2.4. Extend profiles to lower/higher values. In some cases, the
    % sample itself may not have caught the entire profile likelihood. When
    % this happens, the edges of the profile that was made so far
    % (collected in <parprof>) will not be above the chi2 criterion for the
    % single parameters. The strategy is to continue profiling to lower
    % and/or higher values of the parameter.
    %
    % Note: Extending the profile is done in steps of <gridsp>, which is
    % HALF of the distance between the regular profile points. This implies
    % that extension is done on a finer grid. This was a mistake on my
    % side, but I think it is okay to keep it in: extension is being done
    % without the aid of sample points, so maybe a finer grid is a good
    % idea. Will be a bit slow though. Also note that extension to lower
    % values will generally be triggered when the sample is on the lower
    % edge. This is caused by the fact that I shift the profiling grid
    % slightly to higher values to catch the best value exactly. That is
    % also not problematic.
    
    % BLOCK 2.4.1. Start with profiling to lower values of the parameter,
    % if needed. If we are not clearly above the parameter's CI yet, and
    % haven't hit the lower bound ... continu further down until we are, or
    % until we hit a boundary
    flag_disp = 0; % flag to make sure we only display the status on screen once
    while parprof(1,i_p) > bnds_tmp(i_p,1) && parprof(1,end) < mll+chicrit_single+1
        
        if flag_disp == 0
            disp(['  Extending profile for ',names{parnr}, ' to lower parameter values'])
            flag_disp = 1; % set flag to one so we don't continue to print on screen
        end
        
        pmat_tst(ind_fit,1) = parprof(1,1:end-1)'; % use first entry of <parprof> to continue with
        pmat_tst(parnr,1)   = pmat_tst(parnr,1) - gridsp; % continue on the same grid spacing to lower value for this parameter
        pmat_tst(parnr,1)   = max(pmat_tst(parnr,1),bnds_tmp(i_p,1)); % but make sure it is not below lower bound
        p_tst               = pmat_tst(ind_fit_tst,1); % parameter estimates to be fitted for this point
        
        % Since we are extending the profile into a difficult area of
        % parameter space, do two optimisations: rough followed by
        % normal.
        [phat_tst]         = setup_simplex(p_tst,0,pmat_tst); % do a rough optimisation
        [phat_tst,mll_tst] = setup_simplex(phat_tst,1,pmat_tst); % do a normal optimisation
        
        pmat_tst(ind_fit_tst,1) = phat_tst; % replace values in <pmat_tst> with fitted ones
        parprof = [pmat_tst(ind_fit,1)' mll_tst ; parprof]; % and add new entry *on top* of <parprof>
        coll_ok = [coll_ok ; parprof(1,:)]; % also add this additional value to <coll_ok> as they are candidates for extra sampling
        
        if mll_tst < mll % then we even found a better optimum, so collect it
            par_best = parprof(1,:);
            if mll_rem - mll_tst > SETTINGS_OPTIM.real_better && mll - mll_tst > SETTINGS_OPTIM.real_better
                % Only print on screen when the new optimum is really
                % better than the original one, and better than the one
                % that is now the best.
                disp(['  Better optimum found when extending profile for ',names{parnr}, ' down: ',num2str(mll_tst),' (best was ',num2str(mll),')'])
            end
            mll = mll_tst; % update MLL
        end
        
    end
    
    % BLOCK 2.4.2. Do the same thing to higher values of the parameter, if
    % needed. If we are not clearly above the parameter's CI, and haven't
    % hit the upper bound ... continu further up until we are, or until we
    % hit a boundary.
    flag_disp = 0; % flag to make sure we only display the status on screen once
    while parprof(end,i_p) < bnds_tmp(i_p,2) && parprof(end,end) < mll+chicrit_single+1
        
        if flag_disp == 0
            disp(['  Extending profile for ',names{parnr}, ' to higher parameter values'])
            flag_disp = 1; % set flag to one so we don't continue to print on screen
        end
        
        pmat_tst(ind_fit,1) = parprof(end,1:end-1)'; % use last entry of <parprof> to continue with
        pmat_tst(parnr,1)   = pmat_tst(parnr,1) + gridsp; % continue on the same grid spacing to higher value for this parameter
        pmat_tst(parnr,1)   = min(pmat_tst(parnr,1),bnds_tmp(i_p,2)); % make sure it is not above upper bound
        p_tst               = pmat_tst(ind_fit_tst,1); % parameter estimates to be fitted for this point
        
        % Two optimisation rounds.
        [phat_tst]         = setup_simplex(p_tst,0,pmat_tst); % do a rough optimisation
        [phat_tst,mll_tst] = setup_simplex(phat_tst,1,pmat_tst); % do a normal optimisation
        
        pmat_tst(ind_fit_tst,1) = phat_tst; % replace values in <pmat_tst> with fitted ones
        parprof = [parprof;pmat_tst(ind_fit,1)' mll_tst]; % add new entry *on bottom* of <parprof>
        coll_ok = [coll_ok ; parprof(end,:)]; % also add this additional value to <coll_ok> as they are candidates for extra sampling
        
        if mll_tst < mll % then we even found a better optimum, so collect it
            par_best = parprof(end,:);
            if mll_rem - mll_tst > SETTINGS_OPTIM.real_better && mll - mll_tst > SETTINGS_OPTIM.real_better
                % Only print on screen when the new optimum is really
                % better than the original one, and better than the one
                % that is now the best.
                disp(['  Better optimum found when extending profile for ',names{parnr}, ' up: ',num2str(mll_tst),' (best was ',num2str(mll),')'])
            end
            mll = mll_tst; % update MLL
        end
        
    end
    
    coll_prof{i_p,1} = parprof; % collect this <parprof> in the cell array <coll_prof> (for output)
end

%% BLOCK 3. Check if a better optimum was found, and if so, refine it.

if ~isempty(par_best) % if profiling led to a better optimum
    
    disp('  Finished profiling, running a simplex optimisation on the best fit set found ...')
    % Do another optimisation as profiling has detected a better optimum.
    pmat(ind_fit,1) = par_best(1:end-1)'; % update <pmat> with the best value found during profiling
    pfit            = pmat(ind_fit,1); % copy best values to <pfit>, for parameters that are to be fitted
    [phat,mll]      = setup_simplex(pfit,1,pmat); % do a normal optimisation
    % The <phat> and <mll> will be output of this function.

else
    phat = []; % signal <calc_parspace> that we have NOT found a better value
end

%% BLOCK 4. Prune <coll_prof> for later plotting and saving
% The cell array <coll_prof> contains an n-by-6 matrix for each parameter.
% It contains all of the complete parameter sets that make up the profile
% likelihoods. This is handy when we want to resample from a set on the
% profile, to see if hitting parameters bounds for one parameter affects
% the CI of another parameter, and to augment the sample with the profiled
% sets. However, for plotting and saving, we can use a reduced set: only
% the x and y values for the profile (x values are values of the parameter
% of interest, y values are the associated MLL). The first pruning step is
% thus to remove the information on the non-profiled parameters from the
% array.
% 
% The second step of pruning is to remove parts of the profile that we are
% not interested in. Everything beyond the chi2 cut-off for the 95% joint
% CI of the parameters can be removed. Normally, profiling won't go there.
% However, some part of the profile may become irrelevant after we found a
% better optimum, so that needs to be removed. This could leave the profile
% short (few points, which implies low resolution), which will flag
% re-profiling in <calc_parspace>.

flag_short = 0; % if nothing happens, flag is zero
coll_prof_pruned = cell(n_fit,1); % create empty cell array

for i_p = 1:n_fit % run through the fitted parameters
    
    coll_prof_tmp = [coll_prof{i_p}(:,i_p) coll_prof{i_p}(:,end)]; % create a temporary matrix with the x and y values of this profile
    % This is step 1: removing info on other parameters than the profiled one.
    
    ind_outofbnds = coll_prof_tmp(:,1) > bnds_tmp(i_p,2) | coll_prof_tmp(:,1) < bnds_tmp(i_p,1); % indices to parts of the profiles that cross the min-max bounds
    coll_prof_tmp(ind_outofbnds,:) = [];
    
    % In some cases, we have tails on the profiles that extend to very
    % high values of the MLL (because we found a better optimum
    % during/after profiling). These tails can be pruned in step 2.
    if n_fit > 1
        ind_lo = find(coll_prof_tmp(:,2)<mll+chicrit_joint,1,'first'); % index to first entry that is within the cut-off
        ind_hi = find(coll_prof_tmp(:,2)<mll+chicrit_joint,1,'last');  % index to last entry that is within the cut-off
    else
        % Note: an extra 0.5 is needed when we are fitting only 1
        % parameter; otherwise, the profile is cut short (as the joint
        % criterion is the same as the single criterion).
        ind_lo = find(coll_prof_tmp(:,2)<mll+chicrit_joint+0.5,1,'first'); % index to first entry that is within the cut-off
        ind_hi = find(coll_prof_tmp(:,2)<mll+chicrit_joint+0.5,1,'last');  % index to last entry that is within the cut-off
    end
    
    ind_lo = max(1,ind_lo-1); % take one extra set lower if there is one
    ind_hi = min(size(coll_prof_tmp,1),ind_hi+1); % take one extra set higher if there is one
    
    coll_prof_pruned{i_p} = coll_prof_tmp(ind_lo:ind_hi,:); % only keep the ones without any extreme tails
    % Note: this procedure only keeps the sets within the criterion, and
    % adds one lower and higher. This last step is needed for extremely
    % steep profiles.
    
    coll_prof_pruned{i_p}(isinf(coll_prof_pruned{i_p}(:,2)),2) = mll+100; % replace any INF with a very high value
    % Note: this implies that this part of the profile will be plotted as a
    % very steep line. Further, in <calc_parspace>, this will allow
    % interpolation rather than taking the extreme value. This improves the
    % CIs a little bit when a profile is extremely steep.
 
    % If, after pruning, the profile is rather short, flag that, so
    % <calc_parspace> will do an new round of profiling.
    if size(coll_prof_pruned{i_p},1) < 40 % 40 seems good enough
        flag_short = 1;
    end
    
end
