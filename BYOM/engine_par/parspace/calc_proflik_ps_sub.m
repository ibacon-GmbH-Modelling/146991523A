function [parprof,par_best,coll_ok] = calc_proflik_ps_sub(pmat,i_p,coll_all,edges_cloud,bnds_tmp,names,SETTINGS_OPTIM,WRAP)

% This part of calc_proflik_ps was moved to a separate function to
% faciliate working with parfor (parallel for loop).
% 
% Author     : Tjalling Jager
% Date       : September 2020
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

ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
n_fit   = sum(pmat(:,2));       % number of fitted parameters
gridsz  = 50;                   % number of points for the profile over the parameter's range
mll     = coll_all(1,end);      % read the best min log-likelihood from first position of <coll_all>
% NOTEpar: profiling may find a better MLL, which normally would influence
% the behavious of subsequent profiling. However, with parallel processing,
% this cannot be done.
mll_rem = mll; % this is fine as mll_rem is only used for deciding when to print info on screen

chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1);     % criterion for single-parameter CIs
coll_ok     = []; % start with an empty matrix to collect candidate sets that will go through another round of sampling
par_best    = []; % collect best better parameter set when profiling finds a better optimum

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
    
%     waitbar(((i_p-1)*gridsz+i_g)/(n_fit*gridsz),f,['Profiling fitted parameters']) % update progress bar
    
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
    [phat_tst1,mll_tst1] = setup_simplex(p_tst1,0,pmat_tst1,WRAP); % do a rough optimisation
    
    % BLOCK 2.2.3. Generate a profile point from the previous profile point.
    if i_g > 1 % then we have a previous value on the profile line
        % Then also use the previous value to see if propagating that one works better!
        pmat_tst2(ind_fit,1) = parprof(i_g-1,1:end-1)'; % put parameter values for this point from the profile in the temporary <pmat_tst2>
        pmat_tst2(parnr,1)   = par_rng(i_g); % force the value on the current grid point for the parameter to be profiled
        p_tst2               = pmat_tst2(ind_fit_tst,1); % parameter estimates to be fitted for this point
        [phat_tst2,mll_tst2] = setup_simplex(p_tst2,0,pmat_tst2,WRAP); % do a rough optimisation
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
    [phat_tst,mll_tst]      = setup_simplex(phat_tst,1,pmat_tst,WRAP); % do a normal optimisation
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
    [phat_tst2,mll_tst2] = setup_simplex(p_tst2,1,pmat_tst2,WRAP); % do a standard optimisation
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
    [phat_tst]         = setup_simplex(p_tst,0,pmat_tst,WRAP); % do a rough optimisation
    [phat_tst,mll_tst] = setup_simplex(phat_tst,1,pmat_tst,WRAP); % do a normal optimisation
    
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
    [phat_tst]         = setup_simplex(p_tst,0,pmat_tst,WRAP); % do a rough optimisation
    [phat_tst,mll_tst] = setup_simplex(phat_tst,1,pmat_tst,WRAP); % do a normal optimisation
    
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

% coll_prof{i_p,1} = parprof; % collect this <parprof> in the cell array <coll_prof> (for output)