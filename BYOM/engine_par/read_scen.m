function [out_c,kc] = read_scen(type,c,t,glo)

% Usage: [out_c,kc] = read_scen(type,c,t,glo)
%
% This function is used to calculate a concentration, previously prepared
% in glo, or remove scenarios when type is -5. The code is copied from
% make_scen, but glo is no longer global. This saves a lot of speed when
% this function is called from derivatives (and to a lesser effect from
% simplefun)! Furthermore, a number of other changes to increase speed.
% When using an ODE solver, derivatives and this function are called MANY
% times, and therefore, small speed increases there make a huge
% differences. The function make_scen still has the same functionality to
% allow for backwards compatibility.
%
% <type> is used to decide what to return to the calling function: 
%     -1) then we are called from derivatives, and return one c
%     -2) we are called from simplefun and need to return a matrix Tev
%     -3) we are called for fast/slow kinetics and need to return a vector c
% <c>    identifier of the treatment
% <t>    time point or time vector to use for the exposure scenario
% <glo>  the structure with information (normally global); used as input
%        increases speed
%
% Author     : Tjalling Jager
% Date       : September 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM.

int_loc   = glo.int_scen == c; % logical indexing for speed (helps a lot!)
int_coll  = glo.int_coll{int_loc}; % extract the structure with forcing events or spline info
int_type  = glo.int_type(int_loc); % extract the type of interpolation for this scenario

% We are always called to return an exposure profile (or a single
% concentration). I tried a few ways to further increase speed, but at
% least there is no point to if-then-else or move the things called for
% derivatives to the front of the function.

MF = glo.MF; % this is for using read_scen for LPx/EPx calculations as well!

switch int_type
    
    case 1 % forcing with piecewise polynomial
        % NOTE: MF included!
        
        if type == -2
            if numel(t) == 1 && t == -1 % that is the sign that we need to return the time vector of the scenario
                out_c = int_coll.GridVectors{1}; % this is ONLY the time vector for the spline set
            else
                error('A continuous spline cannot work with the analytical solution, so use an ODE solver.')
            end
        else
            out_c = MF * int_coll(t); % use the structure form already provided to calculate exposure concentrations
        end
        
    case 2 % forcing with events (time periods with constant exposure)
        % NOTE: MF included!
        
        switch type
            case -1 % then we are called from derivatives, and return one c
                
                if length(glo.timevar) == 2 && glo.timevar(2) > 0 % then we know which interval to use!
                    out_c = MF * int_coll(glo.timevar(2),2); % simply read value for correct interval
                else
                    out_c = MF * int_coll(find(int_coll(:,1)<=t,1,'last'),2); % and extract correct current forcing
                end
                
            case -2 % we are called from simplefun and need to return a matrix Tev
                out_c = int_coll; % is Tev
                out_c(:,2) = MF * out_c(:,2); % spply MF to the concentrations, not the time vector
                kc    = 0; % no degradation
                
            case -3 % we are called for fast/slow kinetics and need to return a vector c
                out_c = zeros(length(t),1);
                for i = 1:size(int_coll,1)-1
                    out_c(t<int_coll(i+1,1) & t>=int_coll(i,1)) = int_coll(i,2);
                end
                out_c(t>=int_coll(end,1)) = int_coll(end,2); % after last event till end of t
                out_c = MF * out_c; % apply MF to the concentration vector
                
        end
        
    case 3 % forcing as static renewal with first-order disappearance
        % NOTE: MF included!
        
        switch type
            case -1 % then we are called from derivatives, and return one c
                
                kc = int_coll(end,2); % extract the disappearance rate for this scenario
                int_coll = int_coll(1:end-1,:); % remove last line (which contained kc)
                if length(glo.timevar) == 2 && glo.timevar(2) > 0 % then we know which interval to use!
                    ind_i = glo.timevar(2); % index for interval that we need to use
                    c0    = int_coll(ind_i,2); % starting conc of the previous renewal
                    t0    = t-int_coll(ind_i,1); % time since last renewal
                    out_c = MF * c0 * exp(-kc*t0); % calculate concentration for this time point
                else
                    c0 = int_coll(find(int_coll(:,1)<=t,1,'last'),2); % starting conc of the previous renewal
                    t0 = t-int_coll(find(int_coll(:,1)<=t,1,'last'),1); % time since last renewal
                    out_c = MF * c0 * exp(-kc*t0); % calculate concentration for this time point
                end
                
            case -2 % we are called from simplefun and need to return a matrix Tev
                out_c = int_coll(1:end-1,:); % is Tev
                out_c(:,2) = MF * out_c(:,2); % spply MF to the concentrations, not the time vector
                kc    = int_coll(end,2); % degradation rate for this scenario
                
            case -3 % we are called for fast/slow kinetics and need to return a vector c
                kc    = int_coll(end,2); % extract the degradation rate
                out_c = zeros(length(t),1);
                for i = 1:size(int_coll,1)-1 % not last one, as that is just for the kc
                    ind_t = t<int_coll(i+1,1) & t>=int_coll(i,1);
                    out_c(ind_t) = int_coll(i,2) * exp(-kc*(t(ind_t)-int_coll(i,1)));
                end
                ind_t = t>=int_coll(i+1,1); % there might be some more time points in t
                out_c(ind_t) = int_coll(i,2) * exp(-kc*(t(ind_t)-int_coll(i,1)));
                out_c = MF * out_c; % apply MF to the concentration vector
        end
        
    case 4 % linear forcing functions, used with analytical solution
        % NOTE: MF included!
        
        switch type
            case -1 % then we are called from derivatives, and return one c
                
                if length(glo.timevar) == 2 && glo.timevar(2) > 0  % then we know which interval to use!
                    ind_i = glo.timevar(2); % index for interval that we need to use
                    int_coll_tmp = int_coll(ind_i,:); % extract only relevant row from int_coll
                    int_coll_tmp(1,2:3) = MF * int_coll_tmp(1,2:3); % apply MF to the concentrations AND slopes, not the time vector
                    out_c = int_coll_tmp(1,2) + (t-int_coll_tmp(1,1)) * int_coll_tmp(1,3) ;
                else
                    ind_t = find(int_coll(:,1)<=t,1,'last'); % last previous index in exposure scenario
                    int_coll(:,2:3) = MF * int_coll(:,2:3); % apply MF to the concentrations AND slopes, not the time vector
                    out_c = int_coll(ind_t,2) + (t-int_coll(ind_t,1)) * int_coll(ind_t,3) ;
                end
                
            case -2 % we are called from simplefun and need to return a matrix Tev
                
                out_c = int_coll; % is Tev
                out_c(:,2:3) = MF * out_c(:,2:3); % apply MF to the concentrations AND slopes, not the time vector
                kc    = 0; % no degradation
                
            case -3 % we are called for fast/slow kinetics and need to return a vector c
                
                Tev = int_coll;
                if t(end) > Tev(end,1) % do we ask for more points than in Tev?
                    Tev = [Tev; t(end) 0 0]; % add a dummy time point in Tev (from t)
                end
                out_c = zeros(length(t),1);
                for i = 1:size(Tev,1)-1 % run through event episodes
                    ind_t = t<=Tev(i+1,1) & t>=Tev(i,1); % indices for time points in this event
                    % time points on the exact start of each interval may be
                    % calculated twice, which is okay to make sure we also include
                    % first and last time point in there
                    out_c(ind_t) = Tev(i,2) + (t(ind_t)-Tev(i,1)) * Tev(i,3) ;
                end
                % when t is longer than the scenario Tev, the last point in
                % Tev is used (because the scenario definition made sure
                % that there is a last row with slope=0); when t is shorter
                % than Tev, it should work fine (though a bit inefficient)
                out_c = MF * out_c; % apply MF to the concentration vector
                
        end
end


