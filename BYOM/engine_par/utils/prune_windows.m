function Trange_pruned = prune_windows(Trange,Cw,Twin)

% Usage: Trange_pruned = prune_windows(Trange,Cw,Twin)
% 
% Clever trick to reduce the number of windows to calculate from an
% exposure profile. This was worked out by Neil Sherborne: "First, for a
% given exposure profile list all possible time windows (day 0-21, 1-22
% etc.). Then find the largest minimum concentration across all time
% windows. Every window whose maximum is lower than the largest minimum
% cannot be the worst case, since there is a window where exposure is
% greater at all ages."
% 
% This is not advised when there are feedbacks on the 'elimination rate'
% and the pMoA is assimilation, maintenance or growth. A paper on this
% topic is in preparation by Neil. Under these conditions, there may not be
% a unique EPx. Also when modifying the model equations, care is needed. 
% 
% Author     : Tjalling Jager 
% Date       : August 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo

if isfield(glo,'moa') && isfield(glo,'feedb')
    if glo.feedb(2) == 1 && any(glo.moa(1:3)==1)
        warning('off','backtrace')
        warning('Pruning exposure profile windows is NOT recommended for this combination of pMoA and feedbacks!')
        warning('There is a possibility of (limited) higher effects in other windows, so use extra care.')
        disp(' '), warning('on','backtrace')
    end
end

% First, find window with highest minimum concentration.
maxmin_Cw = 0;
for i_T = 1:length(Trange) % run through all time points (start of window)
    
    T = [Trange(i_T);Trange(i_T)+Twin]; % time window of length Twin
    % locate time window in exposure profile
    ind_1  = find(Cw(:,1)>=T(1),1,'first');
    ind_2  = find(Cw(:,1)>=T(2),1,'first');
    
    if ~isempty(ind_2) % then the end of the window is within the total profile
        Cw_tmp = Cw(ind_1:ind_2,:); % extract only the profile that covers the time window
    else % then the end of the window is outside of the profile
        Cw_tmp = Cw(ind_1:end,:); % take the profile as is
    end 
    if Cw_tmp(1,1) > T(1) % if the profile does not start at the exact point where we want to start
        Cw_0   = interp1(Cw(:,1),Cw(:,2),T(1)); % interpolate to the exact point in the profile
        Cw_tmp = cat(1,[T(1) Cw_0],Cw_tmp);     % and add the interpolated point to the profile
    end
    
    if min(Cw_tmp(:,2)) > maxmin_Cw % we have a new maximum minimum
        maxmin_Cw = min(Cw_tmp(:,2));
    end
    
end

% Next, find the windows for which the peak concentration is lower than the
% highest minimum in all windows
Trange_pruned = Trange; % copy the time vector
for i_T = 1:length(Trange) % run through all time points (start of window)
    
    T = [Trange(i_T);Trange(i_T)+Twin]; % time window of length Twin
    % locate time window in exposure profile
    ind_1  = find(Cw(:,1)>=T(1),1,'first');
    ind_2  = find(Cw(:,1)>=T(2),1,'first');
    
    if ~isempty(ind_2) % then the end of the window is within the total profile
        Cw_tmp = Cw(ind_1:ind_2,:); % extract only the profile that covers the time window
    else % then the end of the window is outside of the profile
        Cw_tmp = Cw(ind_1:end,:); % take the profile as is
    end 
    if Cw_tmp(1,1) > T(1) % if the profile does not start at the exact point where we want to start
        Cw_0   = interp1(Cw(:,1),Cw(:,2),T(1)); % interpolate to the exact point in the profile
        Cw_tmp = cat(1,[T(1) Cw_0],Cw_tmp);     % and add the interpolated point to the profile
    end
    
    if max(Cw_tmp(:,2)) < maxmin_Cw % maximum in this window is below highest minimum
        Trange_pruned(i_T) = NaN; % make this window NaN
    end
    
end

Trange_pruned(isnan(Trange_pruned)) = []; % remove the ones we flagged