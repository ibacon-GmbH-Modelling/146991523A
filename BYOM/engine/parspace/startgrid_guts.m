function par = startgrid_guts(par)

% Usage: par = startgrid_guts(par)
%
% This function prepares a parameter structure for GUTS analyses, with
% min-max ranges, as required for the parameter-space explorer
% (<calc_parspace>). As input, it needs to know the parameter structure
% <par>.
%
% Note: unlike openGUTS, this function will calculate different (tighter)
% bounds for some parameters when <kd> is fixed. For any fixed parameter,
% its range will be set to the best value, which automatically affects
% calculation of the other ranges.
% 
% Note: <startgrid_guts> will only use the time vector of the *first* data
% set for its derivation of search ranges.
% 
% Inputs
% <par>  parameter structure as obtained from the script
% 
% Outputs
% <par>  parameter structure with optimised search ranges in min/max columns
%
% Author     : Tjalling Jager
% Date       : March 2021
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <startgrid> code that is
% distributed as part of the Matlab version of openGUTS (see
% http://www.openguts.info). Therefore, this code is distributed under the
% same license as openGUTS (GPLv3). The modifications are only to ensure
% that the code operates in the general BYOM framework.
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

global glo glo2 X0mat DATA % make the data set and initial states global variables

%% BLOCK 1. Preliminary things 
% First some preliminary checks. This section looks whether the parameter
% names <kd> and <mw> are in the parameter vector, which would indicate a
% GUTS analysis. If they are, their locations in the parameter matrix are
% collected in globals.

names = glo2.names; % names of the model parameters (for capturing slow kinetics)
glo.loc_kd = find(strcmp(names,'kd')==1); 
glo.loc_mi = find(strcmp(names,'mw')==1); % location of threshold on which to spot for slow kinetics

if isempty(glo.loc_kd) || isempty(glo.loc_mi) % then we have an incorrect parameter definition (for catching slo kinetics)
    error('Was expecting parameters kd and mw in the parameter structure ...') 
end

%% BLOCK 2. Find min-max time and exposure concentration across all data sets. 

c_max = 0;   % look for maximum peak concentration across all data sets
c_min = inf; % look for the minimum *non-zero* TWA across all data sets
% t = glo.t;
if isfield(glo,'locS')
    t = DATA{1,glo.locS}(2:end,1); % time vector of the FIRST survival data set
elseif isfield(glo,'loc_h')
    t = DATA{1,glo.loc_h}(2:end,1); % time vector of the FIRST survival data set
end
t_max = max(t);      % largest time point in time vector
t_min = min(t(t>0)); % smalles non-zero entry in time vector

for i_d = 1:size(X0mat,2) % run through all scenarios for calibration
    
    c     = X0mat(1,i_d); % extract scenario identifier (or concentration)
    c_v   = c; % if no exposure profile is specified, simply copy c (it is a constant concentration)
    c_twa = c; % if no exposure profile is specified, c is also the TWA
    if isfield(glo,'int_scen') % if it exists: use it to derive current external conc.
        if ismember(c,glo.int_scen) % is c in the scenario range global?
            c_v = make_scen(-3,c,t); % use make_scen again to derive actual exposure concentration vector
            % the -3 lets make_scen know we are calling for a conc. vector
            c_twa = cumtrapz(t,c_v)/t_max; % CHECK CHECK
            c_twa = c_twa(end);
        end
    end
    c_max = max(c_max,max(c_v));
    c_min = min([c_min,min(c_twa(c_twa>0))]); % look for the minimum *non-zero* TWA across all data sets
end

sdit = glo.sel; % selection for SD (1) or IT (2) or mixed (3)
% Make sure that right parameters are fitted for each special case.
if sdit == 1 % for SD ...
    par.Fs = [1 0 1 1 1]; % never fit the threshold spread, set to 1, normal scale
elseif sdit == 2 % for IT ...
    par.bw = [+inf 0 +inf +inf 0]; % never fit the killing rate, set to infinity, normal scale
end

%% BLOCK 3. Find relevant min-max ranges based on the data set.
% Rules of thumb are used to restrict parameter space, based on the data
% set.

if par.hb(2) == 1 % only if <hb> is fitted
    par.hb = [0.01 1 1e-6 0.07 1]; % fit on normal scale
else % otherwise, make the range equal to the best value
    par.hb = [par.hb(1) 0 par.hb(1) par.hb(1) 1];
end

if par.kd(2) == 1 % only if <kd> is fitted
    par.kd(3) = min([log(20)/(5*365)  , -log(1-0.05)/t_max]);
    par.kd(4) = max([log(20)/(0.5/24) , -log(1-0.99)/(0.1*t_min)]);    
    par.kd(5) = 0; % fit on log scale
    par.kd(1) = mean(par.kd(3:4)); % geomean is in statistics toolbox, but that's not really needed anyway
else % otherwise, make the range equal to the best value
    par.kd = [par.kd(1) 0 par.kd(1) par.kd(1) 1];
end

if par.mw(2) == 1 % only if <mw> is fitted
    par.mw(3) = c_min*(1-exp(-par.kd(3)*(4/24)));
    if sdit == 1 % for SD ...
        par.mw(4) = 0.99*c_max; % max as 99% of maximum peak concentration
    else % for IT and mixed ...
        par.mw(4) = 2*c_max; % max as two times of maximum peak concentration
    end
    par.mw(5) = 1; % fit on normal scale
    par.mw(1) = mean(par.mw(3:4));
else % otherwise, make the range equal to the best value
    par.mw = [par.mw(1) 0 par.mw(1) par.mw(1) 1];
end

if par.bw(2) == 1 % only if <bw> is fitted
    par.bw(3) = -log(0.9)/(c_max*t_max); % same as openGUTS
    par.bw(4) = (24^2*0.95)/(par.kd(3)*c_max*exp(-par.kd(3)*t_max*0.5)); 
    par.bw(5) = 0; % fit on log scale
    par.bw(1) = mean(par.bw(3:4));
else % otherwise, make the range equal to the best value
    par.bw = [par.bw(1) 0 par.bw(1) par.bw(1) 1];
end

if par.Fs(2) == 1 % only if <Fs> is fitted
    par.Fs(3:4) = [1 20]; % this range seems good enough
    par.Fs(1) = mean(par.Fs(3:4));
    par.Fs(5) = 0; % fit on log scale
else % otherwise, make the range equal to the best value
    par.Fs = [par.Fs(1) 0 par.Fs(1) par.Fs(1) 1];
end

if isfield(par,'cK') % only if <cK> exists
    if par.cK(2) == 1 % only if <cK> is fitted
        par.cK(3:4) = [0.1*c_min 2*c_max]; % this range seems good enough
        par.cK(1) = mean(par.cK(3:4));
        par.cK(5) = 1; % fit on normal scale
    else % otherwise, make the range equal to the best value
        par.cK = [par.cK(1) 0 par.cK(1) par.cK(1) 1];
    end
end

% For the simplified, reduced GUTS immobility model, we only need an
% additional killing rate for recovery. For now, this gets the same range
% as the regular killing rate.
if isfield(par,'bwr')   
    if par.bwr(2) == 1 % only if <bwr> is fitted
        par.bwr = par.bw;
    else % otherwise, make the range equal to the input value
        par.bwr = [par.bwr(1) 0 par.bwr(1) par.bwr(1) 1];
    end 
end


