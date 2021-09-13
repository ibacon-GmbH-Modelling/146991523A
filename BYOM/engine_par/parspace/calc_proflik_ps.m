function [coll_prof,coll_ok,phat,mll,coll_prof_pruned,flag_short] = calc_proflik_ps(pmat,coll_all,bnds_tmp,~,SETTINGS_OPTIM,WRAP)

% Usage: [coll_prof,coll_ok,phat,mll,coll_prof_pruned,flag_short] = calc_proflik_ps(pmat,coll_all,bnds_tmp,f,SETTINGS_OPTIM,WRAP)
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
% <SETTINGS_OPTIM> settings read from setup_settings by <calc_optim_ps>
% <WRAP>      a wrapper for things that used to be global (e.g., glo and DATA)
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

%% BLOCK 1. Initial things.

global glo2

names = glo2.names; % names of the model parameters

ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
n_fit   = sum(pmat(:,2));       % number of fitted parameters
% gridsz  = 50;                   % number of points for the profile over the parameter's range

% Note: options for rough/detailed simplex optimisations are set in <setup_simplex>

mll     = coll_all(1,end); % read the best min log-likelihood from first position of <coll_all>
% mll_rem = mll; % remember the initial one entered (<mll> may be updated during profiling)
chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1);     % criterion for single-parameter CIs
chicrit_joint  = 0.5 * SETTINGS_OPTIM.crit_table(n_fit,1); % criterion for joint 95% CI or parameters
chicrit_joint  = max(chicrit_joint,chicrit_single+1.5); % makes sure to take a bit more than the single-par cut-off for 1 and 2 fitted parameters

ind_final   = find(coll_all(:,end) < mll + chicrit_joint,1,'last'); % find index to last point still in 95%CI
coll_prof   = cell(n_fit,1); % initialise an empty cell array that will collect the profiles for plotting later
% coll_ok     = []; % start with an empty matrix to collect candidate sets that will go through another round of sampling
% par_best    = []; % collect best better parameter set when profiling finds a better optimum
edges_cloud = [(min(coll_all(1:ind_final,1:end-1)))' (max(coll_all(1:ind_final,1:end-1)))']; % edges of the joint 95% cloud (matrix)

%% BLOCK 2. Do the actual profiling.
% Within BYOMpar, the various parameters are profiled in parallel. This is
% not exactly the same as the BYOM version, while the MLL is not updated
% until all profiling is finished. I don't think that is very problematic.

coll_best = nan(n_fit,n_fit+1);

parfor i_p = 1:n_fit % run through all fitted parameters
    
    [parprof,par_best,coll_ok] = calc_proflik_ps_sub(pmat,i_p,coll_all,edges_cloud,bnds_tmp,names,SETTINGS_OPTIM,WRAP);

    coll_prof{i_p,1} = parprof; % collect this <parprof> in the cell array <coll_prof> (for output)
    if ~isempty(par_best)
        coll_best(i_p,:) = par_best;
    end
    coll_OK{i_p,1}   = coll_ok;
    
end

coll_ok = [];
for i_p = 1:n_fit % run through all fitted parameters
    coll_ok = cat(1,coll_ok,coll_OK{i_p,1}); % put all output for coll_ok together
end    

if all(isnan(coll_best(:,end)))
    par_best = [];
else
    [mll,ind_mll] = min(coll_best(:,end));
    par_best = coll_best(ind_mll,:);
end

%% BLOCK 3. Check if a better optimum was found, and if so, refine it.

if ~isempty(par_best) % if profiling led to a better optimum
    
    disp('  Finished profiling, running a simplex optimisation on the best fit set found ...')
    % Do another optimisation as profiling has detected a better optimum.
    pmat(ind_fit,1) = par_best(1:end-1)'; % update <pmat> with the best value found during profiling
    pfit            = pmat(ind_fit,1); % copy best values to <pfit>, for parameters that are to be fitted
    [phat,mll]      = setup_simplex(pfit,1,pmat,WRAP); % do a normal optimisation
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
