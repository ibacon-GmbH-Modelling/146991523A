%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>, or the explicit
% function(s) in <simplefun.html simplefun.m>. As input, it gets:
%
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
% zero-variate data (used in <byom_bioconc_extra.html
% byom_bioconc_extra.m>).
%
% * Author: Tjalling Jager
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo) 

Xout2    = []; % additional uni-variate output, not used in this example
% zvd      = []; % additional zero-variate output, not used in this example

% Note: if these options are not used, these variables must be defined as
% empty as they are outputs of this function.

% if needed, calculate model values for zero-variate data from parameter
% set; these lines can be removed if no zero-variate data are used
if ~isempty(glo.zvd) % if there are zero-variate data defined (see byom_bioconc_extra)
    zvd       = glo.zvd; % copy zero-variate data structure to zvd
    zvd.ku(3) = par.Piw(1) * par.ke(1); % add model prediction as third value in zvd
else % if there are no zero-variate data defined (as in byom_bioconc_start)
    zvd       = []; % additional zero-variate output, output defined as empty matrix
end

%% Initial settings
% This part extracts optional settings for the ODE solver that can be set
% in the main script (defaults are set in prelim_checks). The useode option
% decides whether to calculate the model results using the ODEs in
% <derivatives.html derivatives.m>, or the analytical solution in
% <simplefun.html simplefun.m>. Using eventson=1 turns on the events
% handling. Also modify the sub-function at the bottom of this function!
% Further in this section, initial values can be determined by a parameter
% (overwrite parts of X0), and zero-variate data can be calculated. See the
% example BYOM files for more information.

useode   = glo.useode; % calculate model using ODE solver (1) or analytical solution (0)
eventson = glo.eventson; % events function on (1) or off (0)
stiff    = glo.stiff; % set to 1 or 2 to use a stiff solver instead of the standard one

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario
% % if needed, extract parameters from par that influence initial states in X0
% Ci0   = par.Ci0(1); % example: parameter for the internal concentration
% X0(2) = Ci0; % put this parameter in the correct location of the initial vector

%% Calculations
% This part calls the ODE solver (or the explicit model in <simplefun.html
% simplefun.m>) to calculate the output (the value of the state variables
% over time). There is generally no need to modify this part. The solver
% ode45 generally works well. For stiff problems, the solver might become
% very slow; you can try ode15s instead.

c  = X0v(1);     % the concentration (or scenario number)
t  = t(:);       % force t to be a row vector (needed when useode=0)

TE = 0; % dummy for time of events

if useode == 1 % use the ODE solver to calculate the solution
    % Note: set options AFTER the 'if useode == 1' as odeset takes
    % considerable calculation time, which is not needed when using the
    % analytical solution. Also note that the global _glo_ is now input to
    % the derivatives function. This increases calculation speed.
    
    % specify options for the ODE solver; feel free to change the
    % tolerances, if you know what you're doing (for some problems, it is
    % better to set them much tighter, e.g., 1e-9)
    options = odeset; % start with default options
    if eventson == 1
        options = odeset(options,'Events',@eventsfun,'RelTol',1e-4,'AbsTol',1e-7); % add an events function and tigher tolerances
    else
        options = odeset(options,'RelTol',1e-4,'AbsTol',1e-7); % only specify tightened tolerances
    end
    % options = odeset(options,'InitialStep',max(t)/1000,'MaxStep',max(t)/100); % specify smaller stepsize
    
    % call the ODE solver (try ode15s for stiff problems, and possibly with for pulsed forcings)
    if isempty(options.Events) % if no events function is specified ...
        switch stiff
            case 0
                [~,Xout] = ode45(@derivatives,t,X0,options,par,c,glo);
            case 1
                [~,Xout] = ode113(@derivatives,t,X0,options,par,c,glo);
            case 2
                [~,Xout] = ode15s(@derivatives,t,X0,options,par,c,glo);
        end
    else % with an events functions ... additional output arguments for events:
        % TE catches the time of an event, YE the states at the event, and IE the number of the event
        switch stiff
            case 0
                [~,Xout,TE,YE,IE] = ode45(@derivatives,t,X0,options,par,c,glo);
            case 1
                [~,Xout,TE,YE,IE] = ode113(@derivatives,t,X0,options,par,c,glo);
            case 2
                [~,Xout,TE,YE,IE] = ode15s(@derivatives,t,X0,options,par,c,glo);
        end
    end
else % alternatively, use an explicit function provided in simplefun
    Xout = simplefun(t,X0,par,c,glo);
end

if isempty(TE) || all(TE == 0) % if there is no event caught
    TE = +inf; % return infinity
end

%% Output mapping
% _Xout_ contains a row for each state variable. It can be mapped to the
% data. If you need to transform the model values to match the data, do it
% here. 

% Xout(:,1) = Xout(:,1).^3; % e.g., do something on first column, like cube it ...

% % To obtain the output of the derivatives at each time point. The values in
% % dXout might be used to replace values in Xout, if the data to be fitted
% % are the changes (rates) instead of the state variable itself.
% % dXout = zeros(size(Xout)); % initialise with zeros
% for i = 1:length(t) % run through all time points
%     dXout(i,:) = derivatives(t(i),Xout(i,:),par,c,glo); 
%     % derivatives for each stage at each time
% end

%% Events function
% Modify this part of the code if _eventson_=1. This subfunction catches
% the 'events': in this case, it looks for the external concentration where
% degradation stops. This function should be adapted to the problem you are
% modelling (this one matches the byom_bioconc_... files). You can catch
% more events by making a vector out of _values_.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

Ct      = par.Ct(1); % threshold external concentration where degradation stops
nevents = 1;         % number of events that we try to catch

value       = zeros(nevents,1); % initialise with zeros
value(1)    = X(1) - Ct;        % thing to follow is external concentration (state 1) minus threshold
isterminal  = zeros(nevents,1); % do NOT stop the solver at an event
direction   = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
