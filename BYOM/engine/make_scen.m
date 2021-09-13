function varargout = make_scen(type,varargin)

% Usage: varargout = make_scen(type,varargin)
%
% This function is used to prepare scenarios with time-varying exposure,
% and to apply it in <derivatives> or <simplefun>. <varargin> can contain
% one or more matrices with scenario definitions (Cw, see examples in the
% GUTS package).
%
% <type> is in various ways in the code below, so watch out. For
% constructing profiles in the global <glo>:
% 
% Type = 1 Creates a piecewise polynomial form of the forcing data points
% supplied. A plot is made to view the splines (cubic hermite spline is
% default). Cw matrices are allowed to contain more than one scenario.
%
% Type = 2 Creates a joint global from of the forcing events supplied
% (varargin contains one or more matrices for Cw). The exposure
% concentration is constant in the period between two events. Cw matrices
% are allowed to contain more than one scenario. This creates block pulses.
%
% Type = 3 Static renewal scenario. Creates a joint global from of the
% forcing events supplied (varargin contains one or more matrices for Cw).
% The exposure concentration decreases exponentially in the period between
% two renewal events. Cw matrices are allowed to contain more than one
% scenario.
%
% Type = 4 Linear interpolation in a time series, to be used with
% analytical solutions (e.g., for scaled damage in reduced GUTS models).
% 
% Negative values for <type> are used to signal specific activities from
% this this function, which are used internally in BYOM (e.g., to clear
% scenarios or to derive exposure concentrations from the globals).
% 
% Note that requesting exposure concentrations from a global scenario
% definition is now also done by <read_scen>. Advantage is that <glo> is no
% longer a global but an input, which makes it much faster when using ODE
% solvers. However, this function also keeps that functionality for now.
% 
% Author     : Tjalling Jager 
% Date       : September 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo % the int_coll and int_scen are now part of the general global structure glo

glo.MF = 1; % this sets the global for the multiplication factor to 1
% Note: this MF is used by read_scen, when calculating LPx/EPx, but it
% needs to be defined in all cases, so it is good to do it here.

meth_inter = 'pchip'; % use pchip when user asks for a continuous spline (continuous avoids hard discontinuities)

varargout = {};
if type < 0 % then we are using this function to calculate a concentration, previously prepared in glo!
    % or remove scenarios when type is -5
    c = varargin{1};
    
    if type == -5 % remove one or more scenarios from the globals
        if c(1) == -1
            if isfield(glo,'int_coll') % first check if there are scenarios defined
                glo = rmfield(glo,'int_coll'); % remove all spline info from glo
                glo = rmfield(glo,'int_scen');
                glo = rmfield(glo,'int_type');
            end
        else
            for i = 1:length(c)
                [~,int_loc] = ismember(c(i),glo.int_scen); % where is c in the scenario range global?
                % int_coll = glo.int_coll{int_loc}; % extract the matrix with forcing events
                % int_type = glo.int_type(int_loc); % extract the type of interpolation for this scenario
                glo.int_coll(int_loc) = [];
                glo.int_type(int_loc) = [];
                glo.int_scen(int_loc) = [];
            end
        end
        return % go back to where you came from and skip rest of function
    end
    
    % [~,int_loc] = ismember(c,glo.int_scen); % where is c in the scenario range global?
    int_loc     = glo.int_scen == c;     % logical indexing for speed (helps a lot!)
    int_coll    = glo.int_coll{int_loc}; % extract the structure with forcing events or spline info
    int_type    = glo.int_type(int_loc); % extract the type of interpolation for this scenario
    t           = varargin{2};
    
    % Code below seems a bit silly. It may be better to apply the LPx or
    % EPx when reading the exposure profile when it is used, rather than
    % creating a new profile every time (and keeping the old one in glo as
    % well). This needs some careful planning though ...
    if type < -5 % then we modify a scenario with an LPx
        switch type
            case -6 % modify scenario c with factor LPx
                LPx = t; % this input argument is now the LPx
                glo.int_coll_rem = glo.int_coll{int_loc};
                switch int_type % which type of time-variable exposure scenario?
                    case 1 % forcing with piecewise polynomial
                        glo.int_coll{int_loc}.Values = glo.int_coll{int_loc}.Values * LPx; % modify it with the factor LPx
                    case 2 % forcing events, constant periods of exposure
                        glo.int_coll{int_loc}(:,2) = glo.int_coll{int_loc}(:,2) * LPx; % modify it with the factor LPx
                    case 3 % forcing events for static renewal
                        glo.int_coll{int_loc}(1:end-1,2) = glo.int_coll{int_loc}(1:end-1,2) * LPx; % modify it with the factor LPx
                    case 4 % linear forcings
                        glo.int_coll{int_loc}(:,2:3) = glo.int_coll{int_loc}(:,2:3) * LPx; % modify conc and slope with the factor LPx
                end
            case -7 % then we repair the change
                glo.int_coll{int_loc} = glo.int_coll_rem; % return the scenario to the original
        end
        return % go back to where you came from and skip rest of function
    end
    
    % Code below is for returning a concentration or a Tev matrix with
    % exposure events. This is now also present in read_scen, which is
    % optimised for speed, so preferred! This is kept here as well for
    % backwards compatibility.
    switch int_type
        
        case 1 % forcing with piecewise polynomial
            
            if type == -2
                if numel(t) == 1 && t == -1 % that is the sign that we need to return the time vector of the scenario
                    varargout{1} = int_coll.GridVectors{1}; % this is the time vector for the spline set
                else
                    error('A continuous spline cannot work with the analytical solution, so use glo.useode=1.')
                end
            else
                c = int_coll(t); % use the structure form already provided
                varargout{1} = c;
            end
            
        case 2 % forcing with events (time periods with constant exposure)
            
            switch type
                case -1 % then we are called from derivatives, and return one c
                    c = int_coll(find(int_coll(:,1)<=t,1,'last'),2); % and extract correct current forcing
                    varargout{1} = c;
                case -2 % we are called from simplefun and need to return a matrix Tev
                    Tev = int_coll;
                    kc  = 0; % no degradation
                    varargout{1} = Tev;
                    varargout{2} = kc;
                case -3 % we are called for fast/slow kinetics and need to return a vector c
                    c = zeros(length(t),1);
                    for i = 1:size(int_coll,1)-1
                        c(t<int_coll(i+1,1) & t>=int_coll(i,1)) = int_coll(i,2);
                    end
                    c(t>=int_coll(end,1)) = int_coll(end,2); % after last event till end of t
                    varargout{1} = c;
            end
            
        case 3 % forcing as static renewal with first-order disappearance
            
            switch type
                case -1 % then we are called from derivatives, and return one c
                    kc = int_coll(end,2); % extract the disappearance rate for this scenario
                    int_coll = int_coll(1:end-1,:); % remove last line (which contained kc)
                    c0 = int_coll(find(int_coll(:,1)<=t,1,'last'),2); % starting conc of the previous renewal
                    t0 = t-int_coll(find(int_coll(:,1)<=t,1,'last'),1); % time since last renewal
                    c = c0 * exp(-kc*t0); % calculate concentration for this time point
                    varargout{1} = c;
                case -2 % we are called from simplefun and need to return a matrix Tev
                    Tev = int_coll(1:end-1,:);
                    kc  = int_coll(end,2); % degradation rate for this scenario
                    varargout{1} = Tev;
                    varargout{2} = kc;
                case -3 % we are called for fast/slow kinetics and need to return a vector c
                    kc = int_coll(end,2); % extract the degradation rate
                    c  = zeros(length(t),1);
                    for i = 1:size(int_coll,1)-1 % not last one, as that is just for the kc
                        ind_t = t<int_coll(i+1,1) & t>=int_coll(i,1);
                        c(ind_t) = int_coll(i,2) * exp(-kc*(t(ind_t)-int_coll(i,1)));
                    end
                    ind_t = t>=int_coll(i+1,1); % there might be some more time points in t
                    c(ind_t) = int_coll(i,2) * exp(-kc*(t(ind_t)-int_coll(i,1)));
                    varargout{1} = c;
            end
            
        case 4 % linear forcing functions, used with analytical solution
            
            switch type
                case -1 % then we are called from derivatives, and return one c
                    
                    ind_t = find(int_coll(:,1)<=t,1,'last'); % last previous index in exposure scenario
                    c = int_coll(ind_t,2) + (t-int_coll(ind_t,1)) * int_coll(ind_t,3) ;
                    varargout{1} = c;
                           
                case -2 % we are called from simplefun and need to return a matrix Tev
                    
                    Tev = int_coll;
                    kc  = 0; % no degradation
                    varargout{1} = Tev;
                    varargout{2} = kc;
                    
                case -3 % we are called for fast/slow kinetics and need to return a vector c
                    
                    Tev = int_coll;
                    if t(end) > Tev(end,1) % do we ask for more points than in Tev?
                        Tev = [Tev; t(end) 0 0]; % add a dummy time point in Tev (from t)
                    end
                    c   = zeros(length(t),1);
                    for i = 1:size(Tev,1)-1 % run through event episodes
                        ind_t = t<=Tev(i+1,1) & t>=Tev(i,1); % indices for time points in this event
                        % time points on the exact start of each interval may be
                        % calculated twice, which is okay to make sure we also include
                        % first and last time point in there
                        c(ind_t) = Tev(i,2) + (t(ind_t)-Tev(i,1)) * Tev(i,3) ;
                    end
                    % when t is longer than the scenario Tev, the last
                    % point in Tev is used (because below the scenario
                    % definition made sure that there is a last row with
                    % slope=0); when t is shorter than Tev, it should work
                    % fine (though a bit inefficient)
                    varargout{1} = c;
                    
            end
    end
    
    return % go back to where you came from and skip rest of function
end

% if they are already defined, extract int_coll and int_scen from glo
if isfield(glo,'int_scen')
    int_coll = glo.int_coll;
    int_scen = glo.int_scen;
    int_type = glo.int_type;
else % otherwise, define as empty
    int_coll = [];
    int_scen = [];
    int_type = [];
end

% have to catch cases were glo may not be completely defined when calling make_pp
if isfield(glo,'basenm') 
    filenm = glo.basenm; % for saving a plot, use the base filename
else
    filenm = 'exposure'; % but if not defined, use a default
end

L_scen = length(int_scen); % how many scenarios we already have

% If varargin contains more than one scenario, we have to do some
% re-arranging (this can be made more efficient in the future, but since
% this is function only called once before an analysis, speed is not an
% issue). The varag2 will contain one scenario in each cell only, even when
% varargin may contain more than one per cell.

% First, we need to make sure that we can deal with scenarios being entered
% as cell arrays (or as matrices, or both in one call!)
vararg_tmp = {};
for i = 1:length(varargin) % run through the input arguments
    if iscell(varargin{i}) % user may have entered a cell array for one or more Cw's
        vararg_tmp = [vararg_tmp varargin{i}]; % add it to the cell array vararg_tmp as cells
    else % otherwise, it is a regular matrix
        vararg_tmp = [vararg_tmp {varargin{i}}]; % so add it to the cell array vararg_tmp as cell
    end
end
varargin = vararg_tmp;

vararg2 = {};
for i = 1:length(varargin)
    if size(varargin{i},2) == 2 % then this cell entry contains just 1 scenario
        vararg2 = [vararg2 varargin{i}]; % add it to the cell array vararg2
    else % more scenarios in this cell entry
        for j = 1:size(varargin{i},2)-1 % run through the scenarios
            vararg2 = [vararg2 varargin{i}(:,[1 j+1])];
        end
    end
end

max_t = 0;
max_c = 0;
max_s = 0;

for i = 1:length(vararg2) % run through all data sets
    int_scen(i+L_scen) = vararg2{i}(1,2); % extract the scenario
    t_scen = vararg2{i}(2:end,1); % extract the time points
    if t_scen(1) ~= 0
        error(['Exposure scenario ',num2str(int_scen(i+L_scen)),': make sure that the scenario definition starts at t=0.'])
        % this might need to be refined when glo.Tinit is used ... however,
        % not starting at zero is probably not common (certainly not in
        % combination with time-varying exposure), and the user should be
        % able to find a workaround by redefining the data set.
    end
end

switch type
    
    case 1 % forcing as a spline, interpolation between all data points
       
        for i = 1:length(vararg2) % run through all data sets
            splpts = vararg2{i}(2:end,1:2); % extract the spline points
            % newer versions of Matlab have different way of making interpolation then pp
            % for old versions of Matlab (before 2011), this will lead to an error
            int_coll{i+L_scen} = griddedInterpolant(splpts(:,1),splpts(:,2),meth_inter,'nearest');
            max_t = max(max_t,max(splpts(:,1))); % look for max of time vector across all sets
            max_c = max(max_c,max(splpts(:,2))); % look for max of conc vector across all sets
            max_s = max(max_s,size(splpts,1)); % look for max number of points across all sets
        end
                
    case 2 % forcing as a series of events with constant concentration
        
        for i = 1:length(vararg2)
            int_coll{i+L_scen} = vararg2{i}(2:end,1:2); % extract the events
            if int_coll{i+L_scen}(1,1) ~= 0
                warning('The first time point of your forcing events in int_coll is NOT zero.')
                warning('This will lead to errors in derivatives, unless your time vector for modelling also does not start at zero.')
            end
            max_t = max(max_t,max(int_coll{i+L_scen}(:,1))); % look for max of time vector across all sets
            max_c = max(max_c,max(int_coll{i+L_scen}(:,2))); % look for max of conc vector across all sets
        end
        
    case 3 % forcing as static-renewal
        
        for i = 1:length(vararg2)
            int_coll{i+L_scen} = vararg2{i}(2:end,1:2); % extract the events
            % note that last row in int_coll is final time and kc only
            
            if int_coll{i+L_scen}(1,1) ~= 0
                warning('The first time point of your forcing events in int_coll is NOT zero.')
                warning('This will lead to errors in derivatives, unless your time vector for modelling also does not start at zero.')
            end
            max_t = max(max_t,max(int_coll{i+L_scen}(:,1))); % look for max of time vector across all sets
            max_c = max(max_c,max(int_coll{i+L_scen}(:,2))); % look for max of conc vector across all sets
        end
        
    case 4 % linear extrapolation in time series, analytical
        
        for i = 1:length(vararg2)
            Cw = vararg2{i}(2:end,1:2); % extract the events            
            Cw(isnan(Cw(:,2)),:) = []; % remove rows with NaNs in the concentration column
            
            % third column will be the linear slope in this interval
            te_coll = [Cw(1:end-1,1:2) diff(Cw(:,2))./diff(Cw(:,1))];
            
            if isempty(te_coll) % this happens when there is only an entry at t=0 (and NaN otherwise)
                te_coll = [Cw 0;vararg2{i}(end,1) Cw(1,2) 0]; % then take constant concentration
                % note that I add a final time point, otherwise it still goes wrong
            end
            te_coll(1+find(diff(te_coll(:,3))==0),:) = []; % prune it: remove values were slope remains the same!
            % Note: the similar algorithm in openGUTS (prepare_data) leaves
            % the last point in. I don't think this is needed here, and I
            % cannot figure out why I did it in openGUTS (as it seems to do
            % nothing there).
            % 
            % Note: pruning needs to be done before removing any INFs, because
            % we allow multiple entries at the same time point. We thus can
            % have two intervals with the same slope that are not continuously
            % increasing (e.g., a sawtooth). This may leave more events in
            % te_coll than strictly necessary.
            
            te_coll(isnan(te_coll(:,3)),:) = []; % simply remove the rows with slope NaN (this allows double time points)
            te_coll(isinf(te_coll(:,3)),:) = []; % simply remove the rows with slope INF (this allows double time points)
            
            if te_coll(end,1) < Cw(end,1) % this will be almost always ... only with a single entry will this be skipped
                te_coll = [te_coll; Cw(end,1:2) 0]; % add a last row to mark end of scenario with slope zero
                % this ensures that extrapolation beyond the scope of the
                % scenario points gives the same value as the last point;
                % the end time of te_coll may be less than the end of the
                % original Cw profile, but that is fine
            end
            
            % TEST TEST extra round of pruning!! This trick can be used in
            % openGUTS as well. The previous pruning round leaves too many
            % events in due to the possibility of NaNs and INFs. The
            % extraneous sets are now removed here.
            %
            % First find end of each interval (extrapolate from start)
            t_end  = te_coll(1:end-1,2) + te_coll(1:end-1,3) .* (te_coll(2:end,1)-te_coll(1:end-1,1));
            t_strt = te_coll(:,2); % collect the start points
            % Find the points where the end of interval i equals start of
            % interval i+1, AND the slopes are equal.
            ind_prn = find(t_end(:) == t_strt(2:end) & te_coll(1:end-1,3) == te_coll(2:end,3));
            te_coll(ind_prn+1,:) = []; % that means we can prune the NEXT interval away
            
            int_coll{i+L_scen} = te_coll;
                 
            % to find the max exposure, we need to consider also the END of
            % each interval (which may not be the start of the next with
            % double timepoints)
            max_t = max(max_t,max(Cw(:,1))); % look for max of time vector across all sets
            max_c = max(max_c,max([t_end;t_strt])); % look for max of conc vector across all sets
        end
        
end

if max_c == 0
    max_c = 1;
end

% put the coll and scen back into the global glo
glo.int_coll = int_coll;
glo.int_scen = int_scen;
glo.int_type = [int_type type*ones(1,length(vararg2))];

if length(unique(int_scen)) ~= length(int_scen) % more splines per scenario requested
    error('You entered multiple data sets for the same scenario. Check your script, and use only one spline per scenario.')
end

if ~isfield(glo,'scen_plot') || glo.scen_plot == 1 % only make a plot if we were not asked not to
        
    t = linspace(0,max_t*1.1,500); % model vector for plotting the scenarios
    
    if max_s > 100 % we have a lot of points to spline ...
        t = linspace(0,max_t*1.05,100*max_s); % model vector for plotting the scenarios
    end
    
    % make one plot with subplots for each state variable
    n = ceil(sqrt(length(vararg2)));
    m = ceil(length(vararg2)/n);
    
    [figh,ft] = make_fig(m,n); % make figure of correct size
    
    for i = 1:length(vararg2) % run through scenarios
        subplot(m,n,i)
        set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        
        hold on
        
        switch type
            case 1
                c = int_coll{i+L_scen}(t);                
                plot(t,c,'k-','LineWidth',2)
                if length(vararg2{i}(2:end,1)) < 50 % for a limited nr of points, plot points as well
                    plot(vararg2{i}(2:end,1),vararg2{i}(2:end,2),'ko','MarkerFaceColor','w','LineWidth',2)
                end
                
                xlabel('time',ft.name,ft.label)
                if isfield(glo,'LabelTable') && ismember(int_scen(i+L_scen),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                    Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == int_scen(i+L_scen)}; % look up the label belonging to the j-th scenario
                    ylabel(['forcing: ',Ltmp],ft.name,ft.label)
                else
                    ylabel(['forcing, scenario ',num2str(int_scen(i+L_scen))],ft.name,ft.label)
                end
                ylim([0 1.02*max_c])
                
            case 2
                % It is easier (and better for error checking) to ask
                % read_scen to reconstruct the exposure scenario, rather
                % than plotting the data points as entered.
                
                te_scen = vararg2{i}(1,2); % extract the scenario
                te_coll = vararg2{i}(2:end,1:2); % extract the events
                
                c = read_scen(-3,te_scen,t,glo); % use read_scen to derive actual exposure concentration vector
                % the -3 lets read_scen know that we need a conc. vector
                
%                 c = zeros(length(t),1);
%                 for j = 1:size(te_coll,1)-1
%                     c(t<te_coll(j+1,1) & t>=te_coll(j,1)) = te_coll(j,2);
%                 end
%                 c(t>=te_coll(end,1)) = te_coll(end,2); % after last event till end of t
                
                plot(t,c,'k-','LineWidth',2)
                plot(te_coll(:,1),te_coll(:,2),'ko','MarkerFaceColor','w','LineWidth',2)
                
                xlabel('time',ft.name,ft.label)
                if isfield(glo,'LabelTable') && ismember(te_scen,glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                    Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == te_scen}; % look up the label belonging to the j-th scenario
                    ylabel(['forcing: ',Ltmp],ft.name,ft.label)
                else
                    ylabel(['forcing, scenario ',num2str(te_scen)],ft.name,ft.label)
                end
                ylim([0 1.02*max_c])
                
            case 3
                % It is easier (and better for error checking) to ask
                % read_scen to reconstruct the exposure scenario, rather
                % than plotting the data points as entered.
                
                te_scen = vararg2{i}(1,2); % extract the scenario
                te_coll = vararg2{i}(2:end,1:2); % extract the events
                
                c = read_scen(-3,te_scen,t,glo); % use read_scen to derive actual exposure concentration vector
                % the -3 lets read_scen know that we need a conc. vector

%                 kc = vararg2{i}(end,2); % extract the degradation rate
%                 c = zeros(length(t),1);
%                 for j = 1:size(te_coll,1)-1 % not last one, as that is just for the kc
%                     ind_t = t<te_coll(j+1,1) & t>=te_coll(j,1);
%                     c(ind_t) = te_coll(j,2) * exp(-kc*(t(ind_t)-te_coll(j,1)));
%                 end
%                 ind_t = t>=te_coll(j+1,1); % there might be some more time points in t for plotting
%                 c(ind_t) = te_coll(j,2) * exp(-kc*(t(ind_t)-te_coll(j,1)));
                
                plot(t,c,'k-','LineWidth',2)
                plot(te_coll(1:end-1,1),te_coll(1:end-1,2),'ko','MarkerFaceColor','w','LineWidth',2)
                
                xlabel('time',ft.name,ft.label)
                if isfield(glo,'LabelTable') && ismember(te_scen,glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                    Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == te_scen}; % look up the label belonging to the j-th scenario
                    ylabel(['forcing: ',Ltmp],ft.name,ft.label)
                else
                    ylabel(['forcing, scenario ',num2str(te_scen)],ft.name,ft.label)
                end
                
                ylim([0 1.02*max_c])
                
            case 4
                % Note: with the changes to this type (allowing double time
                % points and NaNs), it is now easier (and better for error
                % checking) to ask read_scen to reconstruct the exposure
                % scenario, rather than plotting the data points as
                % entered.
                
                te_scen = vararg2{i}(1,2); % extract the scenario
                te_coll = vararg2{i}(2:end,1:2); % extract the events as well
                
                c = read_scen(-3,te_scen,t,glo); % use read_scen to derive actual exposure concentration vector
                % the -3 lets read_scen know that we need a conc. vector
                       
                plot(t,c,'k-','LineWidth',2)
                if size(te_coll,1) < 50 % for a limited nr of points, plot points as well
                    plot(te_coll(:,1),te_coll(:,2),'ko','MarkerFaceColor','w','LineWidth',2)
                end
                
                xlabel('time',ft.name,ft.label)
                if isfield(glo,'LabelTable') && ismember(te_scen,glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                    Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == te_scen}; % look up the label belonging to the j-th scenario
                    ylabel(['forcing: ',Ltmp],ft.name,ft.label)
                else
                    ylabel(['forcing, scenario ',num2str(te_scen)],ft.name,ft.label)
                end
                ylim([0 1.02*max_c])
                
        end
        
        xlim([0 max(t)])
        
    end
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    h_txt = text(0.5, 1,'Exposure scenarios','HorizontalAlignment','center','VerticalAlignment', 'top');
    set(h_txt,ft.name,ft.text); % use standard formatting for this header
    
    % note that plot is not saved is make_pp is called before defining glo.saveplt
    if isfield(glo,'saveplt') && glo.saveplt > 0 % if we want to save the plot
        savenm = ['exp_scen_',filenm];%
        save_plot(figh,savenm,h_txt);
    end
    drawnow

end  



