%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout TE] = call_deri(t,par,X0v,glo)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>. It is specific
% for the model in the directory 'simple_compound', and modified to deal
% with time-varying exposure by splitting up the ODE solving to avoid hard
% switches. 
% 
% As input, it gets:
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_ a vector with initial states and one concentration (scenario number)
% * _glo_ the structure with various types of information (used to be global)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <derivatives.html derivatives.m>. The
% optional output _TE_ is the time at which an event takes place (specified
% using the events function). The events function is set up to catch
% discontinuities. It should be specified according to the problem you are
% simulating. If you want to use parameters that are (or influence) initial
% states, they have to be included in this function. Optional output Xout2
% is for additional uni-variate data (not used here), and zvd is for
% zero-variate data (not used here).
%
% * Author: Tjalling Jager
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, not used in this case
zvd      = []; % additional zero-variate output, not used in this case

%% Initial settings
% This part extracts optional settings for the ODE solver that can be set
% in the main script (defaults are set in prelim_checks). Further in this
% section, initial values can be determined by a parameter (overwrite parts
% of X0), and zero-variate data can be calculated. See the example BYOM
% files for more information.
% 
% This package will always use the ODE solver; therefore the general
% settings in the global glo for the calculation (useode and eventson) are
% removed here. This packacge also always uses the events function, so the
% option eventson is not used. Furthermore, simplefun.m is removed.

stiff      = glo.stiff; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
break_time = glo.break_time; % break time vector up for ODE solver (1) or don't (0)
min_t      = 500; % minimum length of time vector (affects ODE stepsize only, when needed)
names_sep  = glo.names_sep;

if length(stiff) == 1 % second element is used for tolerances
    stiff(2) = 1; % by default: normally tightened tolerances
end

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario
c  = X0v(1);     % the concentration (or scenario number)

% Start with a check on time-varying exposure scenarios. This is a
% time-saver! The isfield and ismember calls take some time that rapidly
% multiplies as derivatives is called many times.
glo.timevar = [0 0]; % flag for time-varying exposure (second element is to tell read_scen which interval we need)
if isfield(glo,'int_scen') && ismember(c,glo.int_scen) % is c in the scenario range of the global?
    glo.timevar = [1 0]; % tell derivatives that we have a time-varying treatment
end

% Deal with fitting multiple data sets with common parameters
if glo.timevar(1) == 1 % this requires scenarios to be used!
    if c >= 100 && ~isempty(names_sep) % then we have more data sets, and more separate parameters per set!
        i_d = floor(c/100); % extract the data set number from the treatment identifier
        for i_sep = 1:length(names_sep) % run through extra parameter names for separate sets
            eval(['par.',names_sep{i_sep},' = par.',[names_sep{i_sep},num2str(i_d)],';']); % copy extra parameter level to par
        end
    end
end

% if needed, extract parameters from par that influence initial states in X0
% start from specified initial size in a model parameter
L0           = par.L0(1); % initial body length (mm) is a parameter
X0(glo.locL) = L0;        % put this estimate in the correct location of the initial vector

%% Calculations
% This part calls the ODE solver to calculate the output (the value of the
% state variables over time). There is generally no need to modify this
% part. The solver ode45 generally works well. For stiff problems, the
% solver might become very slow; you can try ode15s instead.

t     = t(:);   % force t to be a row vector (needed when useode=0)
t_rem = t;      % remember the original time vector (as we will add to it)

% The code below is meant for discontinous time-varying exposure. It will
% also work for constant exposure. Breaking up will not be as effective for
% piecewise polynomials (type 1), and slow, so better turn that off
% manually when you try splining. The code will run each interval between
% two exposure 'events' (in Tev) separately, stop the solver, and restart.
% This ensures that the discontinuities are no problem anymore.

TE  = 0; % dummy for time of events in the events function
Tev = [0 c]; % exposure profile events setting: without anything else, assume it is constant

InitialStep = max(t)/100; % specify initial stepsize
MaxStep     = max(t)/10;  % specify maximum stepsize
% For constant concentrations, and when breaking the time vector, we can
% use a default; Matlab uses as default the length of the time vector
% divided by 10.

if glo.timevar(1) == 1 % time-varying concentrations?
    % if it exists: use it to derive time vector for the events in the exposure scenario    
    Tev = read_scen(-2,c,-1,glo); % use read_scen to derive exposure concentration events
    % the -2 lets read_scen know we need events, the -1 that this is also needed for splines!
    
    min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*2); 
    % For very long exposure profiles (e.g., FOCUS profiles), we now
    % automatically generate a larger time vector (twice the number of the
    % relevant points in the scenario). This is only used to set step size
    % when not using break_time, and for locating maximum length under
    % no-shrinking.
    
    if break_time == 0
        InitialStep = t(end)/(10*min_t); % initial step size
        MaxStep     = t(end)/min_t;      % maximum step size
        % For the ODE solver, when we have a time-varying exposure set
        % here, we base minimum step size on min_t. Small stepsize is a
        % good idea for pulsed exposures; otherwise, stepsize may become so
        % large that a concentration change is missed completely. When we
        % break the time vector, limiting step size is not needed.
    end
end

% This is a means to include a delay caused by the brood pounch in species
% like Daphnia. The repro data are for the appearance of neonates, but egg
% production occurs earlier. This global shifts the model output in this
% function below. This way, the time vector in the data does not need to be
% manipulated, and the model plots show the neonate production as expected.
bp = 0; % by default, no brood-pouch delay
if glo.Tbp > 0 % if there is a need for a brood-pouch delay ... (glo.Tbp must be defined!)
    tbp = t(t>glo.Tbp)-glo.Tbp; % extra times needed to calculate brood-pounch delay
    t   = unique([t;tbp]);      % add the shifted time points to account for brood-pouch delay
    bp  = 1;                    % signal rest of code that we need brood-pouch delay
end

% When an animal cannot shrink in length, we need a long time vector as we
% need to catch the maximum length over time.
if glo.len == 2 && length(t) < min_t % make sure there are at least min_t points
    t = unique([t;(linspace(t(1),t(end),min_t))']);
end

% specify options for the ODE solver
options = odeset; % start with default options for the ODE solver
% This needs further study ... events function removed. Events function
% needs to be considered very carefully for this model (and would only be
% useful for SD).
switch stiff(2)
    case 1 % for ODE15s, slightly tighter tolerances seem to suffice (for ODE113: not tested yet!)
        RelTol  = 1e-4; % relative tolerance (tightened)
        AbsTol  = 1e-7; % absolute tolerance (tightened)
    case 2 % somewhat tighter tolerances ...
        RelTol  = 1e-5; % relative tolerance (tightened)
        AbsTol  = 1e-8; % absolute tolerance (tightened)
    case 3 % for ODE45, very tight tolerances seem to be necessary in some cases
        RelTol  = 1e-9; % relative tolerance (tightened)
        AbsTol  = 1e-9; % absolute tolerance (tightened)
end
options = odeset(options,'RelTol',RelTol,'AbsTol',AbsTol,'Events',@eventsfun,'InitialStep',InitialStep,'MaxStep',MaxStep); % set options
% Note: setting tolerances is pretty tricky. For some cases, tighter
% tolerances are needed but not for others. For ODE45, tighter tolerances
% seem to work well, but not for ODE15s.

T = Tev(:,1); % time vector with events
if T(end) > t(end) % scenario may be longer than t(end)
    T(T>t(end)) = []; % remove all entries that are beyond the last time point
    % this may remove one point too many, but that will be added next
end
if T(end) < t(end) % scenario may (now) be shorter than we need
    T = cat(1,T,t(end)); % then add last point from t
end
% If a lag time is used, it is a good idea to add it as an event as ODE45
% does not like such a switch either.
Tlag = par.Tlag(1);
if Tlag > 0
    T = unique([par.Tlag(1);T]);
end

% BROOD POUCH POINTS WERE HERE, BUT MOVED UP, BEFORE ADDING POINTS IN
% RELATION TO glo.len!

% Breaking up the time vector is faster and more robust when calibrating on
% data from time-varying exposure that contain discontinuities. However,
% when running through very detailed exposure profiles (e.g., those from
% FOCUS), it is probably best to use break_time=0. Testing indicates that
% breaking up will be much, much slower than simply running the ODE solver
% across the entire profile (and produces very similar output).

if break_time == 0 
    
    % simply use the ODE solver for the entire time vector, adding events into the time vector
    t_tmp = unique([T;t;(T(1:end-1)+T(2:end))/2]); % combine T, t, and halfway-T into new time vector
    switch stiff(1)
        case 0
            [tout,Xout,TE,~,~] = ode45(@derivatives,t_tmp,X0,options,par,c,glo);
        case 1
            [tout,Xout,TE,~,~] = ode113(@derivatives,t_tmp,X0,options,par,c,glo);
        case 2
            [tout,Xout,TE,~,~] = ode15s(@derivatives,t_tmp,X0,options,par,c,glo);
    end
    
else
    
    % Transferring the interval to derivatives is a huge time saver! It is
    % not safe to use i as index to Tev since we may add Tlag to T (which
    % is not in Tev); so better calculate it explicitly. Using <find> may
    % be a bit slow, that is why this is done outside of the loop calling
    % the ODE solver.
    ind_Tev = zeros(length(T)-1,1);
    for i = 1:length(T)-1 % run through all intervals between events
        ind_Tev(i) = find(Tev(:,1) <= T(i),1,'last');
    end
    
    % run the ODE solver piece-wise across all exposure events specified
    tout = []; % start with empty output time vector
    Xout = []; % start with empty output state matrix
    for i = 1:length(T)-1 % run through all intervals between events
        t_tmp = [T(i);T(i+1)]; % start with start and end time for this period
        t_tmp = unique([t_tmp;t(t<t_tmp(2) & t>t_tmp(1))]); % and add the time points from t that fit in there
        glo.timevar(2) = ind_Tev(i); % tell derivatives which part of Tev we're in
        
        % use ODE solver to find solution in this interval
        % Note: the switch-case may be compromising calculation speed; however,
        % it needs to be investigated which solver performs best
        switch stiff(1)
            case 0
                [tout_tmp,Xout_tmp,TE,~,~] = ode45(@derivatives,t_tmp,X0,options,par,c,glo);
            case 1
                [tout_tmp,Xout_tmp,TE,~,~] = ode113(@derivatives,t_tmp,X0,options,par,c,glo);
            case 2
                [tout_tmp,Xout_tmp,TE,~,~] = ode15s(@derivatives,t_tmp,X0,options,par,c,glo);
        end
        
        tout = cat(1,tout,tout_tmp); % add output times to previous ones
        Xout = cat(1,Xout,Xout_tmp); % add results (states) to previous ones
        X0   = Xout_tmp(end,:)';     % update X0 for next time window
        % Note: if t_tmp is only 2 elements long, the ODE solver will
        % return all evaluated time points. This is fine as the code below
        % will sieve out the time points to save. Further, each point in T will
        % occur twice, which is also fine. It might be more efficient though to
        % fish out the correct values in the for loop already ...
    end
    
end

if isempty(TE) || all(TE == 0) % if there is no event caught
    TE = +inf; % return infinity
end

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 

%  Since we need to find the maximum of length over time, a large time
%  vector is needed. We cannot use the events function to find the exact
%  point at which length becomes negative, since the events functions only
%  looks at states, and not at the derivatives!
if glo.len == 2 % when animal cannot shrink in length (but does on weight!)
    L     = Xout(:,glo.locL); % take correct state for body lenght
    maxL  = cummax(L); % copy the vector L to maxL and find cumulative maximum
    Xout(:,glo.locL) = maxL; % replace length is Xout with new maximised body length
end

Xout(:,glo.locS) = max(0,Xout(:,glo.locS)); % make sure survival does not get negative
% In some case, it may become just a bit negative, which means zero.

if bp == 1 % if we need a brood-pouch delay ... 
    [~,loct] = ismember(tbp,tout); % find where the extra brood-pouch time points are in the long Xout
    Xbp      = Xout(loct,glo.locR); % only keep the brood-pouch ones we asked for
end
t = t_rem; % return t to the original input vector

% Select the correct time points to return to the calling function
[~,loct] = ismember(t,tout); % find where the requested time points are in the long Xout
Xout     = Xout(loct,:);     % only keep the ones we asked for

if bp == 1 % if we need a brood-pouch delay ...
    [~,loct] = ismember(tbp+glo.Tbp,t); % find where the extra brood-pouch time points SHOULD BE in the long Xout
    Xout(:,glo.locR)    = 0;   % clear the reproduction state variable
    Xout(loct,glo.locR) = Xbp; % put in the brood-pouch ones we asked for
end

% % To obtain the output of the derivatives at each time point. The values in
% % dXout might be used to replace values in Xout, if the data to be fitted
% % are the changes (rates) instead of the state variable itself.
% % dXout = zeros(size(Xout)); % initialise with zeros
% for i = 1:length(t) % run through all time points
%     dXout(i,:) = derivatives(t(i),Xout(i,:),par,c,glo); 
%     % derivatives for each stage at each time
% end
% % 
% % Note: I think this will still work well, even though the ODE solver was
% % run in parts.

%% Events function
% This subfunction catches the 'events': in this case, this one catches
% when the scaled damage exceeds one of the thresholds, and catches the
% switch at puberty.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

% Note: glo needs to be an input to the events function as well, rather
% than a global since the inputs for the events function must match those
% for derivatives.

Lp = par.Lp(1); % length at puberty (mm)
zb = par.zb(1); % effect threshold for the energy budget
zs = par.zs(1); % effect threshold for survival

nevents = 3; % number of events that we try to catch

value    = zeros(nevents,1); % initialise with zeros
value(1) = X(glo.locD) - zb; % follow when scaled damage exceeds the effect threshold for the energy budget
value(2) = X(glo.locD) - zs; % follow when scaled damage exceed the effect threshold for survival
value(3) = X(glo.locL) - Lp; % follow when body length exceeds length at puberty

isterminal = zeros(nevents,1); % do NOT stop the solver at an event
direction  = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
