function [MinColl,MinCI,ind_traits] = calc_epx_window(par_plot,fname_prof,Twin,opt_ecx,opt_conf)

%  Syntax: [MinColl,MinCI,ind_traits] = calc_epx_window(par_plot,fname_prof,Twin,opt_ecx,opt_conf)
%           PART OF ENGINE_PAR
% Calculate EPx/LPx due to an exposure profile, with moving time window.
% This function should work with every TKTD model you throw at it, as long
% as there is at least one state variables indicated with one of the
% dedicated traits (whose position in the state vector is given by:
% <glo.locS>, <glo.locL> or <glo.locR>). The actual calculations of the EPx
% (with CI) is performed by calc_epx.
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

global glo

filenm    = glo.basenm;

X_excl    = opt_ecx.statsup;   % states to suppress from the calculations (e.g., locS)
par_read  = opt_ecx.par_read;  % when set to 1 read parameters from saved set, but do NOT make CIs
batch_epx = opt_ecx.batch_epx; % when set to 1 use batch mode (no output to screen)
notitle   = opt_ecx.notitle;   % set to 1 to suppress titles above plots
start_neg = opt_ecx.start_neg; % set to 1 to start the moving window at minus window width
Feff      = opt_ecx.Feff;      % effect level (>0 en <1), x/100 in LCx (also used here for ECx)
rob_win   = opt_ecx.rob_win;   % set to 1 to use robust EPx calculation for moving time windows, rather than with fzero
rob_rng   = opt_ecx.rob_rng;   % range within which robust EPx is calculated, and number of points
rob_cust  = opt_ecx.rob_cust;  % custom range for ERA purposes
prune_win = opt_ecx.prune_win; % set to 1 to prune the windows to keep the interesting ones

if isempty(opt_conf)
    type_conf = 0; % then we don't need CIs
else
    type_conf = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    type_conf = max(0,type_conf); % if someone uses -1, set it to zero
end

if batch_epx == 0 % don't show warnings if we're in batch mode
    warning('off','backtrace')
    if ~isempty(glo.names_sep)
        warning('You are using separate parameters per data set (with glo.names_sep)')
        warning('Only the FIRST set will be used for EPx (e.g., only f and not f1, f2, etc.')
        disp(' ')
    end
    if rob_win == 0 && isfield(glo,'moa') && isfield(glo,'feedb')
        if glo.feedb(2) == 1 && any(glo.moa(1:3)==1)
            warning('Calculating standard EPx windows is NOT recommended for this combination of pMoA and feedbacks!')
            warning('There is a possibility of (limited) lower EPx values than those reported.')
            warning('It is advisable to use the robust settings with opt_ecx.rob_win = 1.')
            disp(' ')
        end
    end
    warning('on','backtrace')
end

MinCI = []; % make sure it is defined, even when no CIs are made

% see which states are there
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
else 
    rnd = -1;
end

% Load exposure profile from file to use with linear interpolation
Cw     = load(fname_prof);
Tend   = Cw(end,1); % last time point in profile
Trange = 0:floor(Tend); % start the time window every whole day, starting at t=0 of the profile

if isempty(Twin) % if Twin is left empty, we just calculate the entire profile
    Trange = 0;
    Twin = Tend;
end

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

if batch_epx == 0 % don't show waitbar if we're in batch mode
    if rob_win == 0
        if type_conf > 0 % if we make CIs ...
            f = waitbar(0,'Calculating EPx with CIs for moving time window. Please wait.','Name','calc_epx_window.m');
        else
            f = waitbar(0,'Calculating EPx for moving time window. Please wait.','Name','calc_epx_window.m');
        end
    else
        if type_conf > 0 % if we make CIs ...
            f = waitbar(0,'Calculating robust EPx with CIs for moving time window. Please wait.','Name','calc_epx_window.m');
        else
            f = waitbar(0,'Calculating robust EPx for moving time window. Please wait.','Name','calc_epx_window.m');
        end
    end
end
opt_ecx.batch_epx = 1; % set to batch mode so calc_epx provides no output
opt_ecx.par_read  = 0; % when set to 1 read parameters from saved set, but do NOT make CIs
opt_conf.type     = type_conf;
% NOTE: I change the options, but not batch_epx etc.! That way, this
% funtion keeps the setting it was called with, and only calc_epx will go
% into batch mode.

% Xcoll will collect all EPx output, and is initialised with NaNs
for j = 1:length(Feff)
    Xcoll{j} = nan(length(Trange),length(ind_traits));
    if type_conf > 0
        Xlo{j} = nan(length(Trange),length(ind_traits));
        Xhi{j} = nan(length(Trange),length(ind_traits));
    end
end

% Note: the actual calculations of the EPx (with CIs) are performed by
% calc_epx. Therefore, we cannot use parfor again here.

traitcoll = [];
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

    if ~all(Cw_tmp(:,2)==0) % if there is no exposure, there is no effect
        if rob_win == 0 % use standard routine (with fzero)
            [EPx,EPx_lo,EPx_hi,ind_traits_tmp] = calc_epx(par_plot,Cw_tmp,[],opt_ecx,opt_conf,[],rnd);
        else % use robust routine (with fine grid)
            [EPx,EPx_lo,EPx_hi,ind_traits_tmp] = calc_epx_robust(par_plot,Cw_tmp,[],opt_ecx,opt_conf,[],rnd);
        end
        
        traitcoll = unique([traitcoll ind_traits_tmp]); % collect traits that have meaningful output, at some windows
        for i_trt = 1:length(ind_traits_tmp)
            [~,ind_trt] = ismember(ind_traits_tmp(i_trt),ind_traits);
            for j_eff = 1:length(Feff)
                Xcoll{j_eff}(i_T,ind_trt) = EPx{j_eff}(i_trt);
                if type_conf > 0
                    Xlo{j_eff}(i_T,ind_trt) = EPx_lo{j_eff}(i_trt);
                    Xhi{j_eff}(i_T,ind_trt) = EPx_hi{j_eff}(i_trt);
                end
            end
        end
    end
    
end

if batch_epx == 0
    close(f) % close the waiting bar
end

[~,remX] = setdiff(ind_traits,traitcoll); % traits that are not returned by calc_epx since effect is too small
ind_traits(remX) = []; % remove that trait from the trait list
for j_eff = 1:length(Feff)
    Xcoll{j_eff}(:,remX) = []; % remove that trait from the collected values
    if type_conf > 0
        Xlo{j_eff}(:,remX) = []; % remove that trait from the collected values
        Xhi{j_eff}(:,remX) = []; % remove that trait from the collected values
    end
end

if batch_epx == 0 && ~isempty(remX)
    disp(' ')
    for i = 1:length(remX)
        switch remX(i)
            case locS
                disp('For survival, there is insufficient range of effects to calculate all EPx.')
            case locL
                disp('For body length, there is insufficient range of effects to calculate all EPx.')
            case locR
                disp('For reproduction, there is insufficient range of effects to calculate all EPx.')
            case loc_h
                disp('For healthy animals, there is insufficient range of effects to calculate all EPx.')
        end
    end
end

%% Find where the largest effect occurs and how large it is

for i_X = 1:length(ind_traits) % run through traits
    for i_eff = 1:length(Feff) % run through standard multiplication factors
        % find largest effect, and time window where it occurs, for each
        % state at each MF; note that I ignore effects less than 1
        % permille, and round effects more than 99.9% to 100%
        [val_min,ind_min] = min(Xcoll{i_eff}(:,i_X));
        
        MinColl{i_eff}(i_X,:) = [Feff(i_eff) val_min Trange(ind_min)]; % collect MF, minimum and time at which it occurs
        if type_conf > 0
            MinCI{i_eff}(i_X,:)   = [Xlo{i_eff}(ind_min,i_X) Xhi{i_eff}(ind_min,i_X)]; % also collect CI on the minimum
        end
    end
end

if batch_epx == 1  % in batch mode, we can stop here
    return 
end

%% Plot the results
% This makes a multiplot: the different endpoints (traits) are shown in
% rows (EPx plotted), and the different effect levels in columns. Note that
% the x-axis is the time at which the time window starts; the window is
% shifted across the entire profile, until it extends beyond the profile.

n = length(Feff);
m = length(ind_traits);
[figh,ft] = make_fig(m,n,2); % create a figure window of correct size

for i_X = 1:length(ind_traits) % run through traits
    
    for i_eff = 1:length(Feff) % run through effect levels
        
        h_pl = subplot(m,n,(i_X-1)*length(Feff)+i_eff); % make a sub-plot
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
        if i_eff == 1 % only y-label in first column
            if rob_win == 1 % when using the robust routine
                ylab = ['rEPx ' glo.ylab{ind_traits(i_X)}]; % label
            else
                ylab = ['EPx ' glo.ylab{ind_traits(i_X)}]; % label
            end
            ylabel(ylab,ft.name,ft.label)
        else
            set(h_ax,'YTickLabel',[]); % remove tick labels on y-axis
        end
        if i_X == 1 % only title in first row
            title(['Feff = ',num2str(Feff(i_eff))])
        end
        
        if length(Trange) == 1 % have to come up with something special for when there is only 1 time point
            
            % also plot critical lines for safety margins of 10 and 100
            plot([-1 1],[10 10],'k--')
            plot([-1 1],[100 100],'k--')
            plot([Tstrt Tend],[1 1],'k--')
                        
            errorbar(0,Xcoll{i_eff}(:,i_X),Xcoll{i_eff}(:,i_X)-Xlo{i_eff}(:,i_X),Xhi{i_eff}(:,i_X)-Xcoll{i_eff}(:,i_X),'ko-','MarkerFaceColor','k','LineWidth',1)
            ylim([0.1 10000]) % limit y-axis
            xlim([-1 1]) % limit x-axis
            set(gca, 'YScale', 'log')

        else
            
            if type_conf > 0 % then we have CIs to plot
                % Little trick to fill the area between the two curves, to
                % obtain a coloured confidence interval as a band.
                t2  = [Trange';flipud(Trange')]; % make a new time vector that is old one, plus flipped one
                Xin = [Xlo{i_eff}(:,i_X);flipud(Xhi{i_eff}(:,i_X))]; % do the same for the plot line, hi and lo
                fill(t2(~isnan(Xin)),Xin(~isnan(Xin)),'g','LineStyle','none','FaceAlpha',1) % and fill this object
            end
            
            plot(Trange,Xcoll{i_eff}(:,i_X),'k-','LineWidth',1) % plot effects relative to control
            if rob_win == 1 % when using the robust routine, there are hard min-max bounds
                switch rob_rng(4) % when using the robust routine, there are hard min-max bounds
                    case 0 % log scale
                        ylim([0.1 10000]) % limit y-axis
                        set(gca, 'YScale', 'log')
                    case 1 % linear scale
                        ylim([rob_rng(1) rob_rng(2)]) % limit y-axis
                        if rob_rng(1)>0
                            set(gca, 'YScale', 'log')
                        end
                    case 2 % custom scale
                        ylim([rob_cust(1) rob_cust(end)]) % limit y-axis
                        if rob_cust(1)>0
                            set(gca, 'YScale', 'log')
                        end
                end
            else
                ylim([0.1 10000]) % limit y-axis
                set(gca, 'YScale', 'log')
            end
            xlim([Tstrt Tend]) % limit x-axis

            if type_conf > 0 % then we have CIs to plot
                plot(Trange,Xlo{i_eff}(:,i_X),'k:') % plot effects
                plot(Trange,Xhi{i_eff}(:,i_X),'k:') % plot effects
            end
            
            % also plot critical lines for safety margins of 1, 10 and 100
            plot([Tstrt Tend],[100 100],'k--')
            plot([Tstrt Tend],[10 10],'k--')
            plot([Tstrt Tend],[1 1],'k--')
            plot([0 0],[0.1 10000],'k:') % and for start of profile
                        
        end
        
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
    savenm = ['epx_window_plot_',filenm];%
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
        disp('  Lowest EPx confidence intervals');
    case 1
        disp('  Lowest EPx with CIs: Bayesian 95% credible interval');
    case 2
        disp('  Lowest EPx with CIs: 95% pred. likelihood, shooting method');
    case 3
        disp('  Lowest EPx with CIs: 95% pred. likelihood, parspace explorer');
end
disp('================================================================================')
disp('  Lowest EPx (most sensitive time window) for various traits')
disp('  Note that multiple minima may occur; only first one is indicated')
disp('  Also note that most sensitive window may depend on the effect level')
if rob_win == 1
    switch rob_rng(4)
        case 0 % log scale
            disp(['  Robust EPx calculation between: ',num2str(rob_rng(1)),'-',num2str(rob_rng(2)),' (n=',num2str(rob_rng(3)),', log)'])
        case 1 % linear scale
            disp(['  Robust EPx calculation between: ',num2str(rob_rng(1)),'-',num2str(rob_rng(2)),' (n=',num2str(rob_rng(3)),', lin)'])
        case 2 % custom
            disp(['  Robust EPx calculation with custom vector between: ',num2str(rob_cust(1)),'-',num2str(rob_cust(end)),' (n=',num2str(length(rob_cust)),')'])
            % change rob_rng so that 'hitting bounds' below proceeds as expected
            rob_rng(1) = rob_cust(1);
            rob_rng(2) = rob_cust(end);
    end
end
disp('================================================================================')
for i_eff = 1:length(Feff) % run through effect levels
    disp(['EP',num2str(100*Feff(i_eff))])
    for i_X = 1:length(ind_traits)
        switch ind_traits(i_X)
            case locS
                fprintf('  Survival    : ')
            case locL
                fprintf('  Body length : ')
            case locR
                fprintf('  Reproduction: ')
            case loc_h
                fprintf('  Healthy     : ')
        end
                
        % Assume that time window always starts on a full day (the %3.0f
        % only prints whole numbers). Here, I convert the values into
        % strings; this allows me to use < and > for the robust EPx, which
        % has a hard boundary.
        a1 = sprintf('%#.2f',MinColl{i_eff}(i_X,2));
        if type_conf > 0
            a2 = sprintf('%#.2f',MinCI{i_eff}(i_X,1));
            a3 = sprintf('%#.2f',MinCI{i_eff}(i_X,2));
        end
        if rob_win == 1
            if MinColl{i_eff}(i_X,2)==rob_rng(1) 
                a1 = sprintf('<%#.2f',MinColl{i_eff}(i_X,2));
            elseif MinColl{i_eff}(i_X,2)==rob_rng(2)
                a1 = sprintf('>%#.2f',MinColl{i_eff}(i_X,2));
            end
            if type_conf > 0
                if MinCI{i_eff}(i_X,1)==rob_rng(1)
                    a2 = sprintf('<%#.2f',MinCI{i_eff}(i_X,1));
                end
                if MinCI{i_eff}(i_X,2)==rob_rng(2)
                    a3 = sprintf('>%#.2f',MinCI{i_eff}(i_X,2));
                end
            end
        end
        fprintf('%10s at t = %3.0f ',a1,MinColl{i_eff}(i_X,3))
        if type_conf > 0
            fprintf('(%s - %s)',a2,a3)
        end
        fprintf('\n')

    end
end
disp('================================================================================')
diary off

% % minimum window for repro across effect levels
% MinRep = [];
% for i_eff = 1:length(Feff) % run through standard multiplication factors
%     for i_X = 1:length(ind_traits)
%         if ind_traits(i_X) == locR % then we have repro data
%             % collect MF, relative cumul. repro, and start of time window
%             MinRep = cat(1,MinRep,[Feff(i_eff) MinColl{i_eff}(i_X,2) MinColl{i_eff}(i_X,3)]);
%         end
%     end
% end
% % How to find the 'most sensitive time window'? The window with most effect
% % WILL depend on the value of the MF. Furthermore, for different elements
% % of the sample, the most-sensitive window will be somewhere else.
% % Furthermore, it will differ between endpoints.



