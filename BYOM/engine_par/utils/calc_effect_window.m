function [MinColl,MinCI,ind_traits] = calc_effect_window(par_plot,fname_prof,Twin,opt_ecx,opt_conf)
 
%  Syntax: [MinColl,MinCI,ind_traits] = calc_effect_window(par_plot,fname_prof,Twin,opt_ecx,opt_conf)
%           PART OF ENGINE_PAR
% Calculate toxic effects due to an exposure profile, with moving time
% window and different (fixed) multiplication factors. This function should
% work with every TKTD model you throw at it, as long as there is at least
% one state variables indicated with one of the dedicated traits (whose
% position in the state vector is given by: <glo.locS>, <glo.locL> or
% <glo.locR>).
% 
% Some calculation speed can be gained by adding an option that skips
% recalculation of the control response when running through a sample for
% CIs, and when running through different start points for the time window
% (at least, when only tox parameters are fitted).
% 
% <par_plot>   parameter structure for the best-fit curve; if left empty the
%              structure from the saved sample is used
% <fname_prof> filename for the file containing the exposure profile
% <Twin>       length of time window (days)
% <opt_ecx>    options structure for ECx and EPx calculations
% <opt_conf>   options structure for making confidence intervals
% 
% <MinColl> collects, for each state, the MF, minimum of the trait relative
% to the control (thus largest effect), and time at which it occurs.
% <MinCI> collects the CI for the minimum trait level. Structure is
% MinColl{i_MF}(i_X,:). Output of ind_traits is needed to know which state
% is meant with i_X.
% 
% Author     : Tjalling Jager 
% Date       : August 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 X0mat

WRAP.glo  = glo;
WRAP.glo2 = glo2;
% Note that glo will be changed in this function. That does not affect the
% functioning of WRAP, since WRAP is only used for packunpack here.

names     = glo2.names;
filenm    = glo.basenm;
glo_rem   = glo; % remember global glo before we change it

backhaz   = opt_ecx.backhaz;   % parameter name (as string) to set to zero for LCx/LPx calculation to remove background mortality
setzero   = opt_ecx.setzero;   % parameter names (as string array) for extra paramaters to be set to zero
X_excl    = opt_ecx.statsup;   % states to suppress from the calculations (e.g., locS)
MF        = opt_ecx.mf_range;  % range for MFs to make plots
par_read  = opt_ecx.par_read;  % when set to 1 read parameters from saved set, but do NOT make CIs
batch_epx = opt_ecx.batch_epx; % when set to 1 use batch mode (no output to screen)
notitle   = opt_ecx.notitle;   % set to 1 to suppress titles above plots
start_neg = opt_ecx.start_neg; % set to 1 to start the moving window at minus window width
calc_rgr  = opt_ecx.calc_rgr;  % set to 1 to calculate rgr in calc_effect_window (experimental!)
prune_win = opt_ecx.prune_win; % set to 1 to prune the windows to keep the interesting ones

if isempty(opt_conf)
    type_conf = 0; % then we don't need CIs
    use_par_out = 0; % set to 1 to use par_out, as entered in this function, rather than from saved set
else
    type_conf = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    type_conf = max(0,type_conf); % if someone uses -1, set it to zero
    use_par_out = opt_conf.use_par_out; % set to 1 to use par_out, as entered in this function, rather than from saved set
end

if batch_epx == 0 % don't show warnings if we're in batch mode
    warning('off','backtrace')
    if ~isempty(glo.names_sep)
        warning('You are using separate parameters per data set (with glo.names_sep)')
        warning('Only the FIRST set will be used for EPx (e.g., only f and not f1, f2, etc.')
        disp(' ')
    end
    if isfield(glo,'moa') && isfield(glo,'feedb')
        if glo.feedb(2) == 1 && any(glo.moa(1:3)==1)
            warning('Calculating effect windows at fixed MFs is NOT recommended for this combination of pMoA and feedbacks!')
            warning('There is a possibility of (limited) higher effects at lower MFs, so use extra care.')
            disp(' ')
        end
    end
    disp(' '), warning('on','backtrace')
end

MinCI = []; % make sure it is defined, even when no CIs are made

glo.scen_plot = 0; % do not make a plot when calling make_scen

% see which other states are there
locS = [];
locL = [];
locR = [];
if isfield(glo,'locS') && ~ismember(glo.locS,X_excl) % then we have a state of survival
    locS = glo.locS; % collect the location
end
if isfield(glo,'locL') && ~ismember(glo.locL,X_excl) % then we have a state of body length
    locL = glo.locL; % collect the location
end
if isfield(glo,'locR') && ~ismember(glo.locR,X_excl) % then we have a state of reproduction
    locR = glo.locR; % collect the location
end
% And also add the states for the GUTS immobility package. For now, healthy
% only, since for that trait, it is easy to calculate ECx relative to the
% control (for death and immobile, the control is zero). There is a way to
% calculate EPx for death, but that would require summing healthy and
% immobile animals before calculating the effect (or take 1-death), so that
% is a bit more work.
loc_h = [];
if isfield(glo,'loc_h') % then we have a state of healthy
    loc_h = glo.loc_h; % collect the locations
end

ind_traits = [locS locL locR loc_h]; % indices for the traits we want from Xout
% We seem to have no interest for damage ...

% vector with initial values for the states in the simulations
X0mat_tmp    = X0mat(:,1); % take first column for our analysis
X0mat_tmp(1) = 1;          % call our scenario "1"

% If we need CIs, load the best parameter set and the random sample from file
if type_conf > 0 || isempty(par_plot) % also if par_plot is not provided
    [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    if numel(rnd) == 1 % that means that no sample was found
        type_conf = -1; % no need to produce an error, just do analysis without CI
    end
    if isempty(par_plot) % if no par structure was entered in this function ...
        par_plot = par; % simply use the one from the sample file
    end
    if type_conf < 1 || par_read == 1 % then we don't want to make CIs
        type_conf = 0;  % don't make CIs anymore (when triggered by par_read)
        rnd       = []; % make sample empty
    end
end

% if ~isfield(par_plot,'tag_fitted') % apparently, parameters have not been fitted
%     warning('off','backtrace')
%     warning('You did not fit any parameters, so LCx or LPx is based on the values in the initial parameter matrix par.')
%     warning('Any CIs are made from the saved set in the MAT file.')
%     disp(' '), warning('on','backtrace')
% end

% backhaz is by default set to 'hb' in prelim_checks.m for the GUTS package
% identify background hazard and set it to zero
if isempty(backhaz) || ~isfield(par_plot,backhaz) % we need to make background mortality zero
    error('The function calc_effect_window expects a parameter name in opt_ecx.backhaz, which matches a parameter name in your parameter structure that can be set zero to remove background mortality.')
else
    eval(['par_plot.',backhaz,'(1) = 0;']); % set parameter to zero in par_plot
    loc_zero = strcmp(names,backhaz)==1; % where is this parameter in the par structure? (logical indexing)
end
% allow extra parameters to be set to zero, such as initial concentrations
if ~isempty(setzero)
    if ~iscell(setzero) % just to make sure it will be a cell
        setzero = {setzero}; % turn it into a cell array with 1 element
    end
    for i = 1:length(setzero)
        eval(['par_plot.',setzero{i},'(1) = 0;']); % set parameter to zero in par_plot
        loc_zero = loc_zero == 1 | strcmp(names,setzero{i})==1; % add this parameter to loc_zero (logical indexing)
    end
end

% Load exposure profile from file to use with linear interpolation
Cw     = load(fname_prof);
Tend   = Cw(end,1); % last time point in profile
Trange = 0:floor(Tend); % start the time window every whole day, starting at t=0 of the profile

if ~isempty(Twin)
    t = linspace(0,Twin,10); % for a FOCUS profile there should be no
    % need to specify a detailed time vector as the ODE solver will dictate
    % the time step; in any case, if more calculation detail is needed,
    % that should be dealt with in call_deri and not here. Since the
    % calculation in this function does not use any shortcuts, there is no
    % need to have more detail here.
else % if Twin is left empty, we just calculate the entire profile
    t = linspace(0,Tend,10);
    Trange = 0;
    Twin = Tend;
end

% ================== RGR ======================
if calc_rgr == 1 % when calculating intrinsic rate of increase ...
    t  = linspace(0,t(end),200); % take more-detailed time vector as we need to integrate
    t  = t(:); % make sure it is a column vector
    t2 = mean([t(1:end-1) t(2:end)],2); % new averaged time vector
    Th = glo.Tbp; % hatching time is brood-pouch delay (if any)
    rgr_init = 1; % initial guess for the growth rate
end
% ================== RGR ======================

% We may also want to include time windows that start almost at the end of
% the profile, and thus that extend longer than the profile itself. We can
% solve that by adding two time points after the profile that are zero. The
% last point is probably not needed, but it does not hurt either. Note that
% this uses the smallest difference in the time vector for the next point
% after the profile.
Cw = cat(1,Cw,[Cw(end,1)+min(diff(Cw(:,1))) 0;Cw(end,1)+Twin 0]);

% Have the option to start with windows *before* the exposure profile. This
% would be needed if the pulses come very early in the profile (otherwise,
% they would hit juveniles only).
if start_neg == 1 && length(Trange) > 1 % don't do this if Twin was entered as empty
    Cw = cat(1,[-1*Twin 0;-1*min(diff(Cw(:,1))) 0],Cw);
    Trange = Cw(1,1):floor(Tend); % start the time window every whole day, starting at t=0 of the profile
end
Tstrt = Cw(1,1);

% =================== TEST ================================================
if prune_win == 1
    Trange = prune_windows(Trange,Cw,Twin); % prune the range to the interesting windows
end
% =================== TEST ================================================

% figure
% plot(Cw(:,1),Cw(:,2),'k-')

if type_conf > 0 % if we make CIs ...
    
    if use_par_out == 1
        par = par_plot; % then we'll use the input par, rather than the saved one
        % Note: par_plot must already be structured in the main script, such
        % that the fitted parameters match the ones in the saved set, etc.
        % Note: when par_plot is NOT entered, it will have been made equal
        % to par from the saved set.
    end
    
    n_sets   = size(rnd,1); % number of samples from parameter space
    pmat     = packunpack(1,par,0,WRAP); % transform structure *from saved set* into a regular matrix
    % it is better to use the saved par, as there may be differences in the
    % log-setting of parameters between the saved set and the optimised
    % par_out matrix (especially when using the alllog option in
    % calc_slice).
    
    par_comp(par,par_plot,cat(2,backhaz,setzero)) % compare par from input with the one from the MAT file
    ind_fit    = (pmat(:,2)==1); % indices to fitted parameters
    ind_logfit = (pmat(:,5)==0 & pmat(:,2)==1); % indices to pars on log scale that are also fitted!
end

if batch_epx == 0 % don't show waitbar if we're in batch mode
    if type_conf > 0 % if we make CIs ...
        f = waitbar(0,'Calculating moving time window with various MFs and CIs. Please wait.','Name','calc_effect_window.m');
    else
        f = waitbar(0,'Calculating moving time window with various MFs. Please wait.','Name','calc_effect_window.m');
    end
end

% Xcoll will collect output, Xlo/Xhi the CIs, initialised with NaNs
for i = 1:length(MF) % run through standard multiplication factors
    Xcoll{i} = nan(length(Trange),1+length(ind_traits));
    Xlo{i} = nan(length(Trange),length(ind_traits));
    Xhi{i} = nan(length(Trange),length(ind_traits));
end

% Start/check parallel pool
% That is always needed here since I also use the parallel toolbox now for
% the different MFs, even without CIs! This may change if we're going to do
% batch processing.
if glo2.n_cores > 0
    poolobj = gcp('nocreate'); % get info on current pool, but don't create one just yet
    if isempty(poolobj) % if there is no parallel pool ...
        parpool('local',glo2.n_cores) % create a local one with specified number of cores
    end
end

for i_T = 1:length(Trange) % run through all time points (start of window)
    
    if batch_epx == 0
        waitbar(i_T/length(Trange),f) % update the waiting bar
    end
    
    T = [Trange(i_T);Trange(i_T)+Twin]; % time window of length Twin
    % locate time window in exposure profile
    ind_1  = find(Cw(:,1)>=T(1),1,'first');
    ind_2  = find(Cw(:,1)>=T(2),1,'first');
    
    if ~isempty(ind_2) % then the end of the window is within the total profile
        Cw_tmp = Cw(ind_1:ind_2,:); % extract only the profile that covers the time window
    else % then the end of the window is outside of the profile
        Cw_tmp = Cw(ind_1:end,:); % take the profile as is
    end % NOTE: is this still needed with the extension of Cw above?
    if Cw_tmp(1,1) > T(1) % if the profile does not start at the exact point where we want to start
        Cw_0   = interp1(Cw(:,1),Cw(:,2),T(1)); % interpolate to the exact point in the profile
        Cw_tmp = cat(1,[T(1) Cw_0],Cw_tmp);     % and add the interpolated point to the profile
    end
    
    Cw_tmp(:,1) = Cw_tmp(:,1)-T(1); % make time vector for the short profile start at zero again
    Cw_tmp      = [1 1;Cw_tmp];     % add a first row with a scenario identifier
    
    make_scen(-5,-1);    % remove all spline info
    make_scen(4,Cw_tmp); % create the globals to define the forcing function
    
    % First calculate the control response in this time window. It is
    % usually superfluous to do this for each start point of the time
    % window: since there is no exposure in the control, the response will
    % be the same in each window. However, somebody might implement
    % time-varying temperatures or food conditions ...
    Xout   = call_deri(t,par_plot,[0;X0mat_tmp(2:end)],glo); % use call_deri.m to provide the output for the control (concentration zero)
    
    % Note: modifying X0mat_tmp is more efficient than setting glo.MF=0, as
    % in that case, call_deri will still treat it as a time-varying
    % exposure, and run through it in steps.
   
    if calc_rgr == 0 % when NOT calculating intrinsic rate of increase ...
        Xctrl  = Xout(end,ind_traits); % remember the relevant output for the traits (end of time window)
    else % when calculating intrinsic rate of increase ...
        cR = Xout(:,locR); % take cumulative repro from the output
        S  = Xout(:,locS); % take survival from the output
        Xctrl    = rgr_calc(cR,S,t,t2,Th,rgr_init);
        rgr_init = [0 1.2*Xctrl]; % update initial guess (new one will not be far from old one)
    end
    
    % some trickery to get parfor running ...
    glo_tmp   = glo;
    if calc_rgr == 0 % when NOT calculating intrinsic rate of increase ...
        t2        = [];
        Th        = [];
        rgr_init  = [];
    end

    % than go through the multiplication factors
    parfor i_MF = 1:length(MF) % run through standard multiplication factors
        
        Xout = calc_effect_window_sub(MF(i_MF),t,par_plot,X0mat_tmp,glo_tmp);
       
        if calc_rgr == 0 % when NOT calculating intrinsic rate of increase ...
            % collect the effect percentage from the control in the first row
            Xcoll{i_MF}(i_T,:) = [MF(i_MF) Xout(end,ind_traits)./Xctrl];
        else
            cR = Xout(:,locR); % take cumulative repro from the output
            S  = Xout(:,locS); % take survival from the output
            Xcoll{i_MF}(i_T,1:2) = [MF(i_MF) rgr_calc(cR,S,t,t2,Th,rgr_init)/Xctrl];
        end
    end
    
    if type_conf > 0 % if we want CIs ... we'll do it again for each element of the sample
        
        % create a huge matrix to catch effect levels, for every set and every case
        Xcoll_tmp = nan(n_sets,length(MF),length(ind_traits));
        
        % some trickery to get parfor running ... First create a pmat_coll with
        % all sets from the sample!
        pmat_coll = cell(n_sets,1); % cell array to construct pmat for each sample
        for k = 1:n_sets % run through all sets in the sample
            pmat(ind_fit,1) = rnd(k,:); % replace values in pmat with the k-th random sample from the MCMC
            % put parameters that need to be fitted on log scale back on normal
            % scale (as call_deri requires normal scale, in contrast to transfer.m)
            if sum(ind_logfit)>0
                pmat(ind_logfit,1) = 10.^(pmat(ind_logfit,1));
            end
            % Note: pmat is still on normal scale here, but the sample in rnd
            % contains the value on a log scale, if a parameter is fitted on log
            % scale. Call_deri needs normal scale structure par_k.
            pmat(loc_zero,1) = 0; % make selected parameter(s) zero in each set of the sample!
            % this has to be done after the tranformation to normal scale!
            pmat_coll{k} = pmat; % collect it in a huge cell array
        end
        
        % some trickery to get parfor running ...
        glo_tmp   = glo;
        L_i_MF    = length(MF);
        if calc_rgr == 0 % when NOT calculating intrinsic rate of increase ...
            t2        = [];
            Th        = [];
        end
        
        parfor k = 1:n_sets % run through all sets in the sample
            
            par_k = packunpack(2,0,pmat_coll{k},WRAP); % transform parameter matrix into a structure
            
            % Calculate the control response for this parameter set. When only
            % the tox parameters are fitted, this is superfluous. However, we
            % should not make a priori assumptions about how this function will
            % be used! E.g., for GUTS cases, hb may be fitted as well, and we
            % need to have the effect relative to the control for THIS set of
            % the sample.
            Xout   = call_deri(t,par_k,[0;X0mat_tmp(2:end)],glo_tmp); % use call_deri.m to calculate the control response
            
            rgr_init2 = 1; % initial guess for the growth rate (needed here to avoid warnings for parfor)
            if calc_rgr == 0 % when NOT calculating intrinsic rate of increase ...
                Xctrl  = Xout(end,ind_traits); % remember the relevant output for the traits (end of time window)
            else % when calculating intrinsic rate of increase ...
                cR = Xout(:,locR); % take cumulative repro from the output
                S  = Xout(:,locS); % take survival from the output
                Xctrl    = rgr_calc(cR,S,t,t2,Th,rgr_init2);
                rgr_init2 = [0 1.2*Xctrl]; % update initial guess (new one will not be far from old one)
            end
            
            % than go through the multiplication factors
            for i_MF = 1:L_i_MF % run through standard multiplication factors
                
                % use of a sub-function is needed to get parfor to cooperate
                Xout = calc_effect_window_sub(MF(i_MF),t,par_k,X0mat_tmp,glo_tmp);
                
                if calc_rgr == 0 % when NOT calculating intrinsic rate of increase ...
                    X_tmp = Xout(end,ind_traits)./Xctrl; % collect the effect percentage relative to the control
                    
                    % % This is the older version, where min-max were taken here already (no good for Bayes!)
                    % Xlo{i_MF}(i_T,:) = min(Xlo{i_MF}(i_T,:),Xout(end,ind_traits)./Xctrl);
                    % Xhi{i_MF}(i_T,:) = max(Xhi{i_MF}(i_T,:),Xout(end,ind_traits)./Xctrl);
                else
                    cR = Xout(:,locR); % take cumulative repro from the output
                    S  = Xout(:,locS); % take survival from the output
                    X_tmp = rgr_calc(cR,S,t,t2,Th,rgr_init2)/Xctrl; % collect the RGR relative to the control
                    
                    % % This is the older version, where min-max were taken here already (no good for Bayes!)
                    % Xlo{i_MF}(i_T,1) = min(Xlo{i_MF}(i_T,1),rgr_tmp/Xctrl);
                    % Xhi{i_MF}(i_T,1) = max(Xhi{i_MF}(i_T,1),rgr_tmp/Xctrl);
                end
                Xcoll_tmp(k,i_MF,:) = X_tmp; % collect the answer!
            end
        end
        
        % Now find the boundaries of the CIs
        for i_X = 1:length(ind_traits) % run through traits
            for i_MF = 1:length(MF) % run through standard multiplication factors
                if type_conf == 1 % then we're doing Bayes
                    Xlo{i_MF}(i_T,i_X) = prctile(Xcoll_tmp(:,i_MF,i_X),2.5,1);
                    Xhi{i_MF}(i_T,i_X) = prctile(Xcoll_tmp(:,i_MF,i_X),97.5,1);
                else % take min-max
                    Xlo{i_MF}(i_T,i_X) = min(Xcoll_tmp(:,i_MF,i_X),[],1);
                    Xhi{i_MF}(i_T,i_X) = max(Xcoll_tmp(:,i_MF,i_X),[],1);
                end
            end
        end
        
        clear Xcoll_tmp % no need for this large matrix anymore
        
    end
end

if batch_epx == 0
    close(f) % close the waiting bar
end

if calc_rgr == 1 % when calculating intrinsic rate of increase ...
    for i_MF = 1:length(MF) % run through standard multiplication factors
        Xcoll{i_MF}(:,3:end) = []; % remove extraneous columns
        Xlo{i_MF}(:,2:end) = []; % remove extraneous columns
        Xhi{i_MF}(:,2:end) = []; % remove extraneous columns
    end   
    ind_traits = -1;
end

%% Find where the largest effect occurs and how large it is

for i_X = 1:length(ind_traits) % run through traits
    for i_MF = 1:length(MF) % run through standard multiplication factors
        % find largest effect, and time window where it occurs, for each
        % state at each MF; note that I ignore effects less than 1
        % permille, and round effects more than 99.9% to 100%
        [val_min,ind_min] = min(Xcoll{i_MF}(:,i_X+1));
        if val_min > 0.999 % avoid apparent effects due to numerical inaccuracy
            val_min = 1;
            ind_min = 1;
        end
        if val_min < 0.001 % avoid extreme minima which make little sense
            val_min = 0;
        end
        MinColl{i_MF}(i_X,:) = [MF(i_MF) val_min Trange(ind_min)]; % collect MF, minimum and time at which it occurs
        MinCI{i_MF}(i_X,:)   = [Xlo{i_MF}(ind_min,i_X) Xhi{i_MF}(ind_min,i_X)]; % also collect CI on the minimum
    end
end

if batch_epx == 1  % in batch mode, we can stop here
    glo = glo_rem; % return glo to its original value
    return 
end

%% Plot the results
% This makes a multiplot: the different endpoints (traits) are shown in
% rows (plotted as response relative to the control), and the different MFs
% in columns. Note that the x-axis is the time at which the time window
% starts; the window is shifted across the entire profile, until it extends
% beyond the profile.

n = length(MF);
m = length(ind_traits);
[figh,ft] = make_fig(m,n,2); % create a figure window of correct size

for i_X = 1:length(ind_traits) % run through traits
    
    for i_MF = 1:length(MF) % run through standard multiplication factors
        
        h_pl = subplot(m,n,(i_X-1)*length(MF)+i_MF); % make a sub-plot
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        if m>1 && n>1 % only shrink white space when there are more than 1 rows and columns
            p = get(h_pl,'position');
            p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
            set(h_pl, 'position', p);
        end
        
        % create axis labels
        if i_X == length(ind_traits) % only x-label at last row
            xlab = ['start ' glo.xlab]; % label
            xlabel(xlab,ft.name,ft.label)
        else
            set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
        end
        if i_MF == 1 % only y-label in first column
            if ind_traits(i_X) == -1
                ylab = 'relative intrinsic rate';
            else
                ylab = ['rel. ' glo.ylab{ind_traits(i_X)}]; % label
            end
            ylabel(ylab,ft.name,ft.label)
        else
            set(h_ax,'YTickLabel',[]); % remove tick labels on y-axis
        end
        if i_X == 1 % only title in first row
            title(['MF = ',num2str(MF(i_MF))])
        end
        
        if length(Trange) == 1 % have to come up with something special for when there is only 1 time point
            
            % also plot critical lines for 10 and 50% effect
            plot([-1 1],[0.9 0.9],'k--')
            plot([-1 1],[0.5 0.5],'k--')
            
            errorbar(0,Xcoll{i_MF}(:,i_X+1),Xcoll{i_MF}(:,i_X+1)-Xlo{i_MF}(:,i_X),Xhi{i_MF}(:,i_X)-Xcoll{i_MF}(:,i_X+1),'ko-','MarkerFaceColor','k','LineWidth',1)
            ylim([0 1.05]) % limit y-axis
            xlim([-1 1]) % limit x-axis
            
        else
            
            if type_conf > 0 % then we have CIs to plot
                % Little trick to fill the area between the two curves, to
                % obtain a coloured confidence interval as a band.
                t2  = [Trange';flipud(Trange')]; % make a new time vector that is old one, plus flipped one
                Xin = [Xlo{i_MF}(:,i_X);flipud(Xhi{i_MF}(:,i_X))]; % do the same for the plot line, hi and lo
                fill(t2(~isnan(Xin)),Xin(~isnan(Xin)),'g','LineStyle','none','FaceAlpha',1) % and fill this object
            end
            
            plot(Trange,Xcoll{i_MF}(:,i_X+1),'k-','LineWidth',1) % plot effects relative to control
            ylim([0 1.05]) % limit y-axis
            xlim([Tstrt Tend]) % limit x-axis
            
            if type_conf > 0 % then we have CIs to plot
                plot(Trange,Xlo{i_MF}(:,i_X),'k:') % plot effects
                plot(Trange,Xhi{i_MF}(:,i_X),'k:') % plot effects
            end
            
            % also plot critical lines for 10 and 50% effect
            plot([Tstrt Tend],[0.9 0.9],'k--')
            plot([Tstrt Tend],[0.5 0.5],'k--')
            plot([0 0],[0 1.05],'k:') % and for start of profile
        end
        
%         % find largest effect, and time window where it occurs, for each
%         % state at each MF
%         [val_min,ind_min] = min(Xcoll{i_MF}(:,i_X+1));
%         if val_min > 0.999 % avoid apparant minima due to numerical inaccuracy
%             val_min = 1;
%             ind_min = 1;
%         end
%         if val_min < 0.001 % avoid extreme minima which make little sense
%             val_min = 0;
%         end
%         MinColl{i_MF}(i_X,:) = [MF(i_MF) val_min Trange(ind_min)]; % collect MF, minimum and time at which it occurs
%         MinCI{i_MF}(i_X,:)   = [Xlo{i_MF}(ind_min,i_X) Xhi{i_MF}(ind_min,i_X)]; % also collect CI on the minimum
    end
end

if notitle == 0
    switch type_conf
        case 0
            tit = 'No confidence intervals';
        case 1
            tit = 'CIs: Bayesian 95% credible interval';
        case 2
            tit = 'CIs: 95% pred. likelihood, shooting method';
        case 3
            tit = 'CIs: 95% pred. likelihood, parspace explorer';
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    h_txt = text(0.5, 1,tit,'HorizontalAlignment','center','VerticalAlignment', 'top');
end

if glo.saveplt > 0 % if we want to save the plot
    savenm = ['effect_window_plot_',filenm];%
    save_plot(figh,savenm);
end

%% Print results on screen

% We might be called with an entire path to a file, rather than just a
% filename. In that case, only display the filename with the lowest folder
% it is in.
a = strfind(fname_prof,filesep);
if length(a) > 2
    fname_prof = fname_prof(a(end-1):end);
end
    
diary('results.out') % collect output in the diary "results.out"
% disp(' ')
disp('Results from moving time window calculations')
disp(['  Exposure from file: ',fname_prof])
disp(['  Length of window  : ',num2str(Twin)])
switch type_conf
    case 0
        disp('  Largest effect without confidence intervals');
    case 1
        disp('  Largest effect with CIs: Bayesian 95% credible interval');
    case 2
        disp('  Largest effect with CIs: 95% pred. likelihood, shooting method');
    case 3
        disp('  Largest effect with CIs: 95% pred. likelihood, parspace explorer');
end
disp('================================================================================')
disp('  Most sensitive time window for various traits')
disp('  Note that multiple minima may occur; only first one is indicated')
disp('  Also note that most sensitive window may depend on the MF')
disp('  Effect is given as relative value of trait, so small values imply large effect')
disp('================================================================================')
for i_MF = 1:length(MF) % run through standard multiplication factors
    disp(['MF = ',num2str(MF(i_MF))])
    for i_X = 1:length(ind_traits)
        switch ind_traits(i_X)
            case -1
                fprintf('  Intrinisic rate: ')
            case locS
                fprintf('  Survival    : ')
            case locL
                fprintf('  Body length : ')
            case locR
                fprintf('  Reproduction: ')
            case loc_h
                fprintf('  Healthy     : ')
        end
                
        if MinColl{i_MF}(i_X,2) > 0.999
            fprintf('no effects \n')
        else
            % assume that time window always starts on a full day (the %3.0f
            % only prints whole numbers)
            fprintf('%#10.2g at t = %3.0f ',MinColl{i_MF}(i_X,2),MinColl{i_MF}(i_X,3))
            fprintf('(%#1.2g - %#1.2g) \n',MinCI{i_MF}(i_X,1),MinCI{i_MF}(i_X,2))
        end
    end
end
disp('================================================================================')
diary off

% minimum window for repro across MFs
MinRep = [];
for i_MF = 1:length(MF) % run through standard multiplication factors
    for i_X = 1:length(ind_traits)
        if ind_traits(i_X) == locR % then we have repro data
            % collect MF, relative cumul. repro, and start of time window
            MinRep = cat(1,MinRep,[MF(i_MF) MinColl{i_MF}(i_X,2) MinColl{i_MF}(i_X,3)]);
        end
    end
end
% How to find the 'most sensitive time window'? The window with most effect
% WILL depend on the value of the MF. Furthermore, for different elements
% of the sample, the most-sensitive window will be somewhere else.
% Furthermore, it will differ between endpoints.

%% Return the globals to their original states

glo = glo_rem; % return glo to its original value

% ==================================================================================
function rgr = rgr_calc(cR,S,t,t2,Th,rgr_init)

if isempty(S)
    S = ones(size(cR)); % just assume no deaths
end
Rp = diff(cR)./(diff(t)); % estimated repro rate, as difference (on time vector t2)
Sp = mean([S(1:end-1) S(2:end)],2); % mean survival probability in an interval delta t (on time vector t2)
if any(Rp>0) % only if there is repro (otherwise, the initial NaN remains)
    % create an anonymous function with the current values for the traits
    crit = @(rgr) trapz(t2,Rp .* Sp .* exp(-rgr*(t2+Th))) -1; % trapezium rule integration
    % extract 1 as we need to find where the integrated thing is 1
    
    if crit(0) <= 0
        rgr = 0; % don't calculate negative rgr's
    else
        rgr = fzero(crit,rgr_init); % use fzero to find the RGR
    end
end

% ==================================================================================
function Xout = calc_effect_window_sub(MF,t,par_k,X0mat_tmp,glo)

% Small sub-function to get parfor running. This piece of code cannot work
% in the main text due to the way that glo_tmp is changed before using it
% in call_deri. However, in a sub-function, it works fine.

glo.MF = MF; % modify scenario info for this scenario (apply MF)
Xout   = call_deri(t,par_k,X0mat_tmp,glo); % use call_deri.m to provide the output for one scenario
