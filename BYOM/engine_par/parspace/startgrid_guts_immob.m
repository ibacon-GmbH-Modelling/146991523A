function par = startgrid_guts_immob(par)

% Usage: par = startgrid_guts_immob(par)
%
% This function prepares a parameter structure for GUTS-immobility
% analyses, with min-max ranges, as required for the parameter-space
% explorer (<calc_parspace>). As input, it needs to know the parameter
% structure <par>.
% 
% Note: the rules for the ranges are based on those for GUTS. However, it
% is unclear to what extent this also holds for the extended model. This
% needs further testing and tweaking!
%
% Note: unlike openGUTS, this function will calculate different (tighter)
% bounds for some parameters when <kd> is fixed. For any fixed parameter,
% its range will be set to the best value, which automatically affects
% calculation of the other ranges.
% 
% Note: startgrid_guts_immob will only use the time vector of the *first*
% data set for its derivation of search ranges.
% 
% Inputs
% <par>  parameter structure as obtained from the script
% 
% Outputs
% <par>  parameter structure with optimised search ranges in min/max columns
%
% Author     : Tjalling Jager
% Date       : July 2021
% Web support: <http://www.openguts.info>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <startgrid> code that is
% distributed as part of the Matlab prototype of openGUTS (see
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
% names <kd> and <mi*> are in the parameter vector, which they should be
% for a GUTS-immobility analysis. If they are, their locations in the
% parameter matrix are collected in globals.

names = glo2.names; % names of the model parameters (for capturing slow kinetics)
glo.loc_kd = find(strcmp(names,'kr')==1); % let's assume that only damage repair can be very slow

if isfield(par,'mi')
    glo.loc_mi(1) = find(strcmp(names,'mi')==1); % location for single threshold on which to spot for slow kinetics
else
    glo.loc_mi(1) = find(strcmp(names,'mii')==1); % location for immobility threshold on which to spot for slow kinetics
    glo.loc_mi(2) = find(strcmp(names,'mid')==1); % also collect location for death threshold for slow kinetics
end

if isempty(glo.loc_kd) || isempty(glo.loc_mi) % then we have an incorrect parameter definition (for catching slo kinetics)
    error('Was expecting parameters kr and mi* in the parameter structure ...') 
end

if isfield(glo,'lim_k')
    lim_k = glo.lim_k;
else
    lim_k = 0;
end
if isfield(glo,'mw_log')
    mw_log = glo.mw_log;
else
    mw_log = 0;
end

%% BLOCK 2. Find min-max time and exposure concentration across all data sets. 

c_max = 0;   % look for maximum peak concentration across all data sets
c_min = inf; % look for the minimum *non-zero* TWA across all data sets
% t = glo.t;
t = DATA{1,glo.loc_h}(2:end,1); % time vector of the FIRST healthy-animals data set
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
% Removed code: this is already taken care of in the main script (unused
% parameters are simply not defined in par).

%% BLOCK 3. Find relevant min-max ranges based on the data set.
% Rules of thumb are used to restrict parameter space, based on the data
% set.

if par.hb(2) == 1 % only if <hb> is fitted
    par.hb = [0.01 1 1e-6 0.07 1]; % fit on normal scale
else % otherwise, make the range equal to the best value
    par.hb = [par.hb(1) 0 par.hb(1) par.hb(1) 1];
end

% NOTE: might need to limit ke as calculation may get very slow otherwise
if par.ke(2) == 1 && glo.fastslow(1) == 0 % only if <ke> is fitted, and under 'normal kinetics'
    par.ke(3) = min([log(20)/(5*365)  , -log(1-0.05)/t_max]);
    par.ke(4) = max([log(20)/(0.5/24) , -log(1-0.99)/(0.1*t_min)]);  
    par.ke(4) = min(par.ke(4),200); % but don't make it too large
    
    if lim_k == 1 % then limit ke for faster analysis
        par.ke([3 4]) = [0.01 10];
    end
    
    par.ke(5) = 0; % fit on log scale
    par.ke(1) = mean(par.ke(3:4)); % geomean is in statistics toolbox, but that's not really needed anyway
else % otherwise, make the range equal to the best value
    par.ke = [par.ke(1) 0 par.ke(1) par.ke(1) 1];
end

% NOTE: might need to limit kr as calculation may get very slow otherwise
if par.kr(2) == 1 && glo.fastslow(2) == 0 % only if <kr> is fitted, and under 'normal kinetics'
    par.kr(3) = min([log(20)/(5*365)  , -log(1-0.05)/t_max]);
    par.kr(4) = max([log(20)/(0.5/24) , -log(1-0.99)/(0.1*t_min)]); 
    par.kr(4) = min(par.kr(4),200); % but don't make it too large
    
    if lim_k == 1 % then limit kr for faster analysis
        par.kr([3 4]) = [0.01 10];
    end
    
    par.kr(5) = 0; % fit on log scale
    par.kr(1) = mean(par.kr(3:4)); % geomean is in statistics toolbox, but that's not really needed anyway
else % otherwise, make the range equal to the best value
    par.kr = [par.kr(1) 0 par.kr(1) par.kr(1) 1];
end

% Estimate mi as in openGUTS, using kr for now!! (that will generally be
% lowest)
mi_est(3) = c_min*(1-exp(-par.kr(3)*(4/24)));
if sdit == 1 % for SD ...
    mi_est(4) = 0.99*c_max; % max as 99% of maximum peak concentration
else % for IT and mixed ...
    mi_est(4) = 2*c_max; % max as two times of maximum peak concentration
end
if mw_log == 0
    mi_est(5) = 1; % fit on normal scale
else
    mi_est(5) = 0; % fit on log scale
end
mi_est(1) = mean(mi_est(3:4));
mi_est(2) = 1;

% And then look which threshold parameters need to be filled with mw_est
if isfield(par,'mi') 
    if par.mi(2) == 1 % only if <mi> is fitted
        par.mi = mi_est;
    else % otherwise, make the range equal to the input value
        par.mi = [par.mi(1) 0 par.mi(1) par.mi(1) 1];
    end 
end
if isfield(par,'mii') 
    if par.mii(2) == 1 % only if <mwi> is fitted
        par.mii = mi_est;
    else % otherwise, make the range equal to the input value
        par.mii = [par.mii(1) 0 par.mii(1) par.mii(1) 1];
    end 
end
if isfield(par,'mid') 
    if par.mid(2) == 1 % only if <mid> is fitted
        par.mid = mi_est;
    else % otherwise, make the range equal to the input value
        par.mid = [par.mid(1) 0 par.mid(1) par.mid(1) 1];
    end 
end

% Estimate bi* as in openGUTS, also using kr for now!!
bi_est(3) = -log(0.9)/(c_max*t_max); % same as openGUTS
bi_est(4) = (24^2*0.95)/(par.kr(3)*c_max*exp(-par.kr(3)*t_max*0.5));
bi_est(5) = 0; % fit on log scale
bi_est(1) = mean(bi_est(3:4));
bi_est(2) = 1;

% And then look which effect-strength parameters need to be filled with bi_est
if isfield(par,'bi')    
    if par.bi(2) == 1 % only if <bi> is fitted
        par.bi = bi_est;
    else % otherwise, make the range equal to the input value
        par.bi = [par.bi(1) 0 par.bi(1) par.bi(1) 1];
    end 
end
if isfield(par,'bii')    
    if par.bii(2) == 1 % only if <bii> is fitted
        par.bii = bi_est;
    else % otherwise, make the range equal to the input value
        par.bii = [par.bii(1) 0 par.bii(1) par.bii(1) 1];
    end 
end
if isfield(par,'bid')    
    if par.bid(2) == 1 % only if <bid> is fitted
        par.bid = bi_est;
    else % otherwise, make the range equal to the input value
        par.bid = [par.bid(1) 0 par.bid(1) par.bid(1) 1];
    end 
end
if isfield(par,'bir')    
    if par.bir(2) == 1 % only if <bir> is fitted
        par.bir = bi_est;
    else % otherwise, make the range equal to the input value
        par.bir = [par.bir(1) 0 par.bir(1) par.bir(1) 1];
    end 
end

% estimate Fs* as in openGUTS
Fs_est(3:4) = [1 50]; % this range seems good enough (note that upper value is larger than in openGUTS!)
Fs_est(1) = mean(Fs_est(3:4));
Fs_est(5) = 0; % fit on log scale
Fs_est(2) = 1;

% And then look which spread parameters need to be filled with Fs_est
if isfield(par,'Fs')
    if par.Fs(2) == 1 % only if <Fs> is fitted
        par.Fs = Fs_est;
    else % otherwise, make the range equal to the input value
        par.Fs = [par.Fs(1) 0 par.Fs(1) par.Fs(1) 1];
    end
end
if isfield(par,'Fsi')    
    if par.Fsi(2) == 1 % only if <Fsi> is fitted
        par.Fsi = Fs_est;
    else % otherwise, make the range equal to the input value
        par.Fsi = [par.Fsi(1) 0 par.Fsi(1) par.Fsi(1) 1];
    end 
end
if isfield(par,'Fsd')    
    if par.Fsd(2) == 1 % only if <Fsd> is fitted
        par.Fsd = Fs_est;
    else % otherwise, make the range equal to the input value
        par.Fsd = [par.Fsd(1) 0 par.Fsd(1) par.Fsd(1) 1];
    end 
end

