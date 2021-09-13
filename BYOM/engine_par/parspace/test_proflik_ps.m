function [flag_profile,coll_ok] = test_proflik_ps(pmat,coll_prof,coll_all,SETTINGS_OPTIM)

% Usage: [flag_profile,coll_ok] = test_proflik_ps(pmat,coll_prof,coll_all,SETTINGS_OPTIM)
% 
% Testing of the profile likelihood, using the sample of sets. This
% function looks for two situations: 
% 1) Points where the MLL of the sample is considerably higher (less good)
%    than the MLL of the profile. This indicates that the sampling has been
%    poor, so we need to resample later.
% 2) Points where the MLL of the sample is considerably lower (better fit)
%    than the profile. This indicates that the profiling was not good
%    enough and may need to be redone.
% 
% Inputs
% <pmat>      parameter matrix
% <coll_prof> profile likelihoods with all parameters for each point (cell array)
% <coll_all>  the entire sample (matrix)
% <SETTINGS_OPTIM> = settings read from setup_settings by <calc_optim_ps>
% 
% Output
% <flag_profile> two flags for checking whether reprofiling is needed (at
%    this moment, a vector with the number of points flagged and the 
%    maximum MLL distance that was found for these points)
% <coll_ok>      matrix to collect candidate parameter sets for resampling
% 
% Author     : Tjalling Jager
% Date       : May 2020
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <test_proflik> code that is
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

ind_fit        = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
coll_ok        = []; % initialise empty matrix to collect candidates for resampling
flag_profile   = [0 0]; % initialise vector for flags with zeros
chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1); % criterion for single-parameter CIs
mll            = coll_all(1,end); % maximum likelihood value so far

%% BLOCK 2. Perform the checks, running through all fitted parameters.
% We run through slices of the parameter space again. The profile
% likelihoods are regularly spaced. As we step through the profiled points,
% we consider the slice from parameter space that is halfway between the
% previous point and the next in the profile.

for i_p = 1:length(ind_fit) % run through all fitted parameters
    
    % BLOCK 2.1. Initial things
    parprof = coll_prof{i_p,1}; % read profile information for this parameter from the ones entered (matrix with 6 columns)
    
    % There can be (extreme) cases where a profile for one parameter is
    % entirely on the wrong place (when profiling has led to a better
    % optimum, but extension of a certain parameter was not triggered).
    % This section checks whether the current profile in <parprof> captures
    % the relevant part of the current sample in <coll_all>. If not, a
    % reprofiling is triggered.
    coll_tst = coll_all(coll_all(:,end) < mll+chicrit_single+1,:); % make a test matrix with the part that is most relevant (not too far from inner rim)
    if min(coll_tst(:,i_p)) < min(parprof(:,i_p)) || max(coll_tst(:,i_p)) > max(parprof(:,i_p)) 
        flag_profile(1) = flag_profile(1) + 1; % count another problem
        flag_profile(2) = +inf; % make sure to flag this as a big problem!
    end
    % Note: there are some cases where this part will be triggered again
    % and again. To solve this, more detail in the profile would be needed.
    % For now, this is 'solved' by the fact that only 2 reprofiles are
    % allowed before stopping the extra sampling.
    
    parprof = parprof(parprof(:,end) < mll+chicrit_single+1,:); % only keep the part that is most relevant (not too far from inner rim)
    % This avoids refinement in areas that we're not interested in.
            
    for i_g = 2:size(parprof,1)-1 % run through grid, but skip first and last (when it hits an edge, strange things can happen)
        
        % BLOCK 2.2. Find best value from <coll_all> in this slice. This
        % parameter set will be compared to the profile likelihood. Note
        % that <coll_all> is sorted, so the first within this range is
        % automatically the best.
        %
        % Note: the spacing between the profile points is not necessarily
        % the same everywhere. Profiles may hit a boundary (and are then
        % placed on top of that boundary) or may be extended at half the
        % default grid spacing. Therefore, it is better to take half the
        % distance between the adjacent points. The <gridsp> is thus a
        % two-element vector with the half distance to the previous point
        % and to the next one. Since grid spacing is smaller when the
        % profile has been extended, it will probably trigger more gaps
        % than is strictly needed, and thus too much resampling. Only
        % downside is a larger sample and longer calculation times, but
        % that should be acceptable (since this will be a tough parameter
        % space anyway).
        
        gridsp = 0.5 * diff(parprof([i_g-1 i_g i_g+1],i_p)); % vector with half distance to previous and next point!
        ind_tst = find(coll_all(:,i_p)>(parprof(i_g,i_p)-gridsp(1)) & coll_all(:,i_p)<(parprof(i_g,i_p)+gridsp(2)),1,'first');
        if isempty(ind_tst) % it can be empty if there are no sample points in this sub-range (which would be very odd)
            [~,ind_tst] = min(abs(coll_all(:,i_p) - parprof(i_g,i_p))); % just take the closest point available
            mll_compare = +inf; % and set MLL here to INF (so it will be flagged for a gap)
        else
            mll_compare = coll_all(ind_tst,end); % this is the best MLL from the sample in this slice
        end
        mll_prof = parprof(i_g,end); % this is the MLL from the profile

        if mll_prof < mll_compare % then the profile is better than the sample (as it should be)
            
            if mll_compare - mll_prof > SETTINGS_OPTIM.gap_extra % then the gap is too wide for comfort
                coll_ok = [coll_ok ; coll_all(ind_tst,:) ; parprof(i_g,:)]; % put the sample point and the profile point in <coll_ok> for refinement
            end
            
        elseif mll_compare < mll_prof % then the sample is better than the profile (which implies a need for new profiling)
            % The profile is on a limited number of grid points on the
            % parameter axis, while the sample is anywhere in the slice.
            % This can be done complex, but it is probably enough to only
            % look at the minimum of the current profile point, previous
            % one, and next one. If the sample is below the lowest of the
            % three, it deserves to be counted as a true deviation.
            %
            % Note: one could decide to also put these problem areas into
            % <coll_ok>. I had that in a previous version as a test, but
            % decided not to do it anymore (can't remember why, though).
            
            min_MLL = min(parprof([i_g-1 i_g i_g+1],end)); % this is the smallest of previous, current and next MLL of the profile
            
            if mll_compare < min_MLL % if the sample is better than the lower edge of the profile in this slice ...
                 flag_profile(1) = flag_profile(1) + 1; % count another problem
                 flag_profile(2) = max(flag_profile(2),min_MLL-mll_compare); % and collect the degree of difference if it is bigger than the previous
            end
        end
    end
end
