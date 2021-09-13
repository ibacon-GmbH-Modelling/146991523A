function [EPx,EPx_lo,EPx_hi,ind_traits] = calc_epx(par_plot,fname_prof,Twin,opt_ecx,opt_conf,opt_tktd,varargin)

% Usage: [EPx,EPx_lo,EPx_hi,ind_traits] = calc_epx(par_plot,fname_prof,Twin,opt_ecx,opt_conf,opt_tktd,varargin)
% 
% Calculate EPx for all available traits with confidence intervals, for a
% given time window. This function should work with every TKTD model you
% throw at it, as long as there is at least one state variables indicated
% with one of the dedicated traits: <glo.locS>, <glo.locL>, <glo.locR>
% (more may be added in the future).
% 
% Here, fzero is used to calculate the EPx, with rough calculations
% providing the interval to search on. This is dangerous since there are
% cases where EPx is not unique. This happens in DEBtox analyses, for
% specific feedback configurations and pMoAs. The calc_epx_robust is better
% suited to catch the correct lowest EPx in those cases, but it is much
% slower.
% 
% Some calculation speed can be gained by adding an option that skips
% recalculation of the control response when running through a sample for
% CIs. At least, when only the tox parameters are fitted!
% 
% <par_plot>   parameter structure for the best-fit curve; if left empty the
%            structure from the saved sample is used
% <fname_prof> filename for the file containing the exposure profile
%              OR exposure profile as two-column matrix!
% <Twin>       time window as two-element vector (empty to use full profile)
% <opt_ecx>    options structure for ECx and EPx calculations
% <opt_conf>   options structure for making confidence intervals
% <opt_tktd>   options structure for plotting results (response at the EPx)
% 
% <EPx> collects, for each effect level and each state, the EPx. <EPx_lo>
% and <EPx_hi> collect the CI for EPx. Structure is EPx{i_F}(i_X). Output
% of ind_traits is needed to know which state is meant with i_X.
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

% When called from calc_epx_window, rnd will be input. This is handy as
% this function will then be called many times (for each window), and it is
% not efficient to load the same rnd from file every time.
if ~isempty(varargin)
    rnd = varargin{1};
else
    rnd = -1;
end

names     = glo2.names;
% filenm    = glo.basenm;

X0mat_rem = X0mat; % remember X0mat before we change it
glo_rem   = glo; % remember glo before we change it

backhaz   = opt_ecx.backhaz; % parameter name (as string) to set to zero for LCx/LPx calculation to remove background mortality
setzero   = opt_ecx.setzero; % parameter names (as string array) for extra paramaters to be set to zero
Feff      = opt_ecx.Feff;    % effect level (>0 en <1), x/100 in LCx (also used here for ECx)
ECx_plot  = opt_ecx.plot;    % set to 0 to NOT make a plot of effects vs time at MFs
X_excl    = opt_ecx.statsup; % states to suppress from the calculations (e.g., locS)
par_read  = opt_ecx.par_read; % when set to 1 read parameters from saved set, but do NOT make CIs
batch_epx = opt_ecx.batch_epx; % when set to 1 use batch mode (no output to screen)

if isempty(opt_conf)
    type_conf   = 0; % then we don't need CIs
    use_par_out = 0; % by default, set to off
else
    type_conf   = opt_conf.type; % use values from slice sampler (1), likelihood region (2) to make intervals
    type_conf   = max(0,type_conf); % if someone uses -1, set it to zero
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
            warning('Calculating standard EPx windows is NOT recommended for this combination of pMoA and feedbacks!')
            warning('There is a possibility of (limited) lower EPx values than those reported.')
            warning('It is advisable to use the robust settings with opt_ecx.rob_win = 1.')
            disp(' ')
        end
    end
    warning('on','backtrace')
end

glo.scen_plot = 0; % do not make a plot when calling make_scen

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
% only, since for that trait, it is easy to calculate EPx relative to the
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

make_scen(-5,-1); % remove all spline info for exposure profiles (just to be on the safe side and perhaps to save some memory)

if numel(rnd)==1 && rnd == -1 % only do this when rnd is not in the input
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
else
    par = par_plot;
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
    error('The function calc_ecx expects a parameter name in opt_ecx.backhaz, which matches a parameter name in your parameter structure that can be set zero to remove background mortality.')
else
    eval(['par_plot.',backhaz,'(1) = 0;']); % set parameter to zero in par_plot
    loc_zero = strcmp(names,backhaz)==1; % where is this parameter in the par structure? (logical indexing)
    % this will be used to make every set in the sample hb=0, so no need to modify par
end
% allow extra parameters to be set to zero, such as initial concentrations
if ~isempty(setzero)
    if ~iscell(setzero) % just to make sure it will be a cell
        setzero = {setzero}; % turn it into a cell array with one element
    end
    for i = 1:length(setzero)
        eval(['par_plot.',setzero{i},'(1) = 0;']); % set parameter to zero in par_plot
        loc_zero = loc_zero == 1 | strcmp(names,setzero{i})==1; % add this parameter to loc_zero  (logical indexing)
    end
end

% Take exposure profile from file or input to use with linear interpolation
if ischar(fname_prof) % if it is a character string ...
    Cw = load(fname_prof); % load from defined file
else
    Cw = fname_prof; % then the exposure profile is entered into this function
end

if isempty(Twin) % if it is empty, just take the full profile
    Cw_tmp = [1 1;Cw]; % add a first row with a scenario identifier
else
    
    % We may also want to include time windows that start almost at the end of
    % the profile, and thus that extend longer than the profile itself. We can
    % solve that by adding two time points after the profile that are zero. The
    % last point is probably not needed, but it does not hurt either.
    Cw = cat(1,Cw,[Cw(end,1)+min(diff(Cw(:,1))) 0;Cw(end,1)+diff(Twin) 0]);
    
    % locate time window in exposure profile
    ind_1 = find(Cw(:,1)>=Twin(1),1,'first');
    ind_2 = find(Cw(:,1)>=Twin(2),1,'first');
    
    if ~isempty(ind_2) % then the end of the window is within the total profile
        Cw_tmp = Cw(ind_1:ind_2,:); % extract only the profile that covers the time window
    else % then the end of the window is outside of the profile
        Cw_tmp = Cw(ind_1:end,:); % take the profile as is
    end % NOTE: is this still needed with the extension of Cw above?
    if Cw_tmp(1,1) > Twin(1) % if the profile does not start at the exact point where we want to start
        Cw_0   = interp1(Cw(:,1),Cw(:,2),Twin(1)); % interpolate to the exact point in the profile
        Cw_tmp = cat(1,[Twin(1) Cw_0],Cw_tmp);     % and add the interpolated point to the profile
    end
    
    Cw_tmp(:,1) = Cw_tmp(:,1)-Twin(1); % make time vector for the short profile start at zero again
    Cw_tmp      = [1 1;Cw_tmp];        % add a first row with a scenario identifier
end

make_scen(4,Cw_tmp); % create the globals to define the forcing function (always linear interpolation)
t = linspace(0,Cw_tmp(end,1),100); % for a FOCUS profile there should be no
% need to specify a detailed time vector as the ODE solver will dictate the
% time step; in any case, if more calculation detail is needed, that should
% be dealt with in call_deri and not here. Since the calculation in this
% function does not use any shortcuts, there is no need to have more detail
% here.
% 
% Note: however, for truly pulsed exposure, make sure to break up the time
% vector in call_deri.

%% Rough exploration of the concentration range
% To find out if there is an ECx,t, and where it approximately is. This
% should give us good starting ranges for all traits and all time points.

if batch_epx == 0
    f = waitbar(0,'Calculating EPx. Please wait.','Name','calc_epx.m');
end

% first calculate the control response in this time window
Xout   = call_deri(t,par_plot,[0;X0mat_tmp(2:end)],glo); % use call_deri.m to provide the output for the control (concentration zero)
Xctrl  = Xout(end,ind_traits); % remember the relevant output for the traits
% Note: modifying X0mat_tmp is more efficient than setting glo.MF=0, as in
% that case, call_deri will still treat it as a time-varying exposure, and
% run through it in steps.

MF        = 1; % start with multplication factor 1
glo.MF    = MF; % modify global
Xout      = call_deri(t,par_plot,X0mat_tmp,glo); % use call_deri.m to provide the output for one scenario
Xout_coll = [MF Xout(end,ind_traits) ./ Xctrl]; % remember the relative output for the traits
% Xout_coll a two columns at the start for multiplication factor

while ~all(min(Xout_coll(:,2:end),[],1)<1-max(Feff)) && MF < 1e6 % stop increasing MF until there is large enough effect for all traits
    MF = MF * 10;
    glo.MF = MF; % modify global
    Xout   = call_deri(t,par_plot,X0mat_tmp,glo); % use call_deri.m to provide the output for one scenario
    Xout   = [MF Xout(end,ind_traits) ./ Xctrl]; % remember the relative output for the traits
    Xout_coll = cat(1,Xout_coll,Xout);
end

MF = 1; % start again from multiplication factor 1
while ~all(max(Xout_coll(:,2:end),[],1)>1-min(Feff)) && MF > 1e-7 % stop decreasing MF until there is small enough effect for all traits
    MF = MF / 10;
    glo.MF = MF; % modify global
    Xout   = call_deri(t,par_plot,X0mat_tmp,glo); % use call_deri.m to provide the output for one scenario
    Xout   = [MF Xout(end,ind_traits) ./ Xctrl]; % remember the relative output for the traits
    Xout_coll = cat(1,Xout,Xout_coll);
end

if MF == 1e-7
    error('It appears that there are effects at MFs much lower than 1; either the risk is very high or something has gone wrong!')
end

% see if this ranges catches all effect levels
if batch_epx == 0
    disp(' ')
end
remX = [];
for i_X = 1:length(ind_traits)
    if min(Xout_coll(:,1+i_X)) > 1-max(Feff) || max(Xout_coll(:,1+i_X)) < 1-min(Feff)
        if batch_epx == 0
            switch ind_traits(i_X) 
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
        if all(Xout_coll(:,1+i_X) > 1-min(Feff))
            % now only remove a state when, at no point at all, is there
            % enough effect for the smallest effect level. In other cases,
            % there may be at least enough to calculate some effect levels.
            remX = cat(2,remX,i_X); % remember that trait for removal
            if batch_epx == 0
                disp('   State is removed from output.')
            end
        end
    end
end

ind_traits(remX)    = []; % remove that trait from the trait list
Xout_coll(:,1+remX) = []; % remove that trait from the collected values
Xctrl(remX)         = []; % remove that trait from the control values

%% Calculate EPx exactly with fzero
% Here, fzero is used to calculate the EPx, with the rough calculations
% providing the interval to search on. This is dangerous since there are
% cases where EPx is not unique. This happens in DEBtox analyses, for
% specific feedback configurations and pMoAs. The calc_epx_robust is better
% suited to catch the correct lowest EPx in those cases, but it is much
% slower.
% 
% This is now also done with the parallel toolbox, running the
% traits/effect levels in parallel. For small numbers of traits/effects
% this will not be faster ... but it is good to speed up the slowest
% analyses.

% Start/check parallel pool
if glo2.n_cores > 0
    poolobj = gcp('nocreate'); % get info on current pool, but don't create one just yet
    if isempty(poolobj) % if there is no parallel pool ...
        parpool('local',glo2.n_cores) % create a local one with specified number of cores
    end
end

% Use parfor already here (one loop with all combinations of traits
% and effect levels). Xout_coll is not large, so don't worry about
% broadcasting it to the parallel toolbox.
XF       = allcomb(1:length(ind_traits),1:length(Feff)); % combine indices for traits and effects
EPx_coll = nan(size(XF,1),1); % initialise matrix to catch EPx
parfor i = 1:size(XF,1)
    i_X = XF(i,1); % read correct index for trait
    i_F = XF(i,2); % read correct index for effects
    
    ind1     = find(Xout_coll(:,i_X+1)<1-Feff(i_F),1,'first');
    EP_range = Xout_coll([ind1-1 ind1],1); % range where EPx,t is located
    if numel(EP_range) == 2
        % use fzero to zero in on the exact value
        EPx_coll(i) = fzero(@calc_epx_sub,EP_range,[],t,par_plot,Feff(i_F),X0mat_tmp,Xctrl(i_X),ind_traits(i_X),glo); % find the EPx,t
    else
        EPx_coll(i) = NaN; % then a proper range was not found
    end
end

% now modify EPx_coll back to the original format of EPx
for i = 1:size(XF,1)
    i_X = XF(i,1);
    i_F = XF(i,2);
    EPx{i_F}(i_X) = EPx_coll(i);
end
clear EPx_coll

% % This is the analysis without the parallel toolbox
% for i_X = 1:length(ind_traits) % run through traits
%     for i_F = 1:length(Feff) % run through effect levels
%         
%         if batch_epx == 0
%             waitbar(((i_X-1)*length(Feff)+i_F-0.5)/(length(ind_traits)*length(Feff)),f) % make a nice waiting bar
%         end
%         
%         ind1    = find(Xout_coll(:,i_X+1)<1-Feff(i_F),1,'first');
%         EP_range = Xout_coll([ind1-1 ind1],1); % range where EPx,t is located
%         if numel(EP_range) == 2
%             % use fzero to zero in on the exact value
%             EPx{i_F}(i_X) = fzero(@calc_epx_sub,EP_range,[],t,par_plot,Feff(i_F),X0mat_tmp,Xctrl(i_X),ind_traits(i_X),glo); % find the EPx,t
%         else
%             EPx{i_F}(i_X) = NaN; % then a proper range was not found
%         end
%         
%     end
% end

if batch_epx == 0
    close(f) % close the waiting bar
end

%% Display results without CI on screen

if type_conf > 0 && batch_epx == 0 % only do this when we'll go into CI calculation next
    
    % disp(' ')
    disp('Results for EPx without CIs (they are calculated next)')
    disp('================================================================================')
    for i_F = 1:length(Feff) % run through effect levels
        disp(['EP',num2str(100*Feff(i_F))])
        
        for i_X = 1:length(ind_traits) % run through traits
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
            fprintf('%#10.2f \n',EPx{i_F}(i_X))
            
        end
        disp('================================================================================')
    end
end

%% Calculate confidence intervals on the EPx

% initialise matrices to catch highest and lowest results from
% sample, per state and per time point
for i_F = 1:length(Feff)
    EPx_lo{i_F} = nan(1,length(ind_traits));
    EPx_hi{i_F} = nan(1,length(ind_traits));
end

if type_conf > 0 % if we make CIs ...
       
    if batch_epx == 0
        disp('Calculation of CIs on EPx with parallel toolbox (no progress will be shown)')
    end
    
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

    % create a huge matrix to catch EPx for every set and every case
    EPx_coll = nan(n_sets,length(Feff),length(ind_traits));
    
    % some trickery to get parfor running ...
    glo_tmp   = glo;
    L_i_F     = length(Feff);
    L_i_X     = length(ind_traits);
    
    parfor k = 1:n_sets % run through all sets in the sample
        
        par_k = packunpack(2,0,pmat_coll{k},WRAP); % transform parameter matrix into a structure
        
        % Calculate the control response for this parameter set. When only
        % the tox parameters are fitted, this is superfluous. However, we
        % should not make a priori assumptions about how this function will
        % be used! E.g., for GUTS cases, hb may be fitted as well, and we
        % need to have the effect relative to the control for THIS set of
        % the sample.
        Xout   = call_deri(t,par_k,[0;X0mat_tmp(2:end)],glo_tmp); % use call_deri.m to provide the output for one scenario
        Xctrl  = Xout(end,ind_traits); % remember the relevant output for the traits
        
        for i_X = 1:L_i_X % run through traits
            for i_F = 1:L_i_F % run through effect levels
                if ~isnan(EPx{i_F}(i_X))
                    % use fzero to zero in on the exact value
                    EPx_tmp = fzero(@calc_epx_sub,EPx{i_F}(i_X),[],t,par_k,Feff(i_F),X0mat_tmp,Xctrl(i_X),ind_traits(i_X),glo_tmp); % find the EPx,t
                    EPx_coll(k,i_F,i_X) = EPx_tmp; % collect the answer!
                end
            end
        end
        
    end
    clear pmat_coll % clear this large variable as it is no longer needed

    % Now find the boundaries of the CIs
    for i_X = 1:length(ind_traits) % run through traits
        for i_F = 1:length(Feff) % run through effect levels
            if type_conf == 1 % then we're doing Bayes
                EPx_lo{i_F}(i_X) = prctile(EPx_coll(:,i_F,i_X),2.5,1);
                EPx_hi{i_F}(i_X) = prctile(EPx_coll(:,i_F,i_X),97.5,1);
            else % take min-max
                EPx_lo{i_F}(i_X) = min(EPx_coll(:,i_F,i_X),[],1);
                EPx_hi{i_F}(i_X) = max(EPx_coll(:,i_F,i_X),[],1);
            end
        end
    end
    clear EPx_coll % no need for this large matrix anymore

end

%% See if we can return now

if batch_epx == 1
    % return the globals to their original states
    X0mat = X0mat_rem;
    glo   = glo_rem;
    return
end

%% Display results on screen

diary(glo.diary) % collect output in the diary "results.out"
% disp(' ')
disp('Results from EPx calculations')
if ischar(fname_prof)
    % We might be called with an entire path to a file, rather than just a
    % filename. In that case, only display the filename with the lowest folder
    % it is in.
    a = strfind(fname_prof,filesep);
    if length(a) > 2
        fname_prof = fname_prof(a(end-1):end);
    end
    disp(['  Exposure from file: ',fname_prof])
else
    disp('  Exposure profile entered directly')
end
switch type_conf
    case 0
        disp('  EPx without confidence intervals');
    case 1
        disp('  EPx with CIs: Bayesian 95% credible interval');
    case 2
        disp('  EPx with CIs: 95% pred. likelihood, shooting method');
    case 3
        disp('  EPx with CIs: 95% pred. likelihood, parspace explorer');
end
disp('================================================================================')
for i_F = 1:length(Feff) % run through effect levels
    disp(['EP',num2str(100*Feff(i_F))])
    
    for i_X = 1:length(ind_traits) % run through traits
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
        
        fprintf('%#10.2f ',EPx{i_F}(i_X))
        if type_conf > 0
            fprintf('(%#8.2f - %#8.2f) ',EPx_lo{i_F}(i_X),EPx_hi{i_F}(i_X))
        end
        fprintf('\n')
        
    end
    disp('================================================================================')
end

diary off

%% Make a plot for the various MFs that are the EPx
% Use plot_tktd to make a series of plots.

if ~isempty(opt_tktd) && ECx_plot ~= 0

    % Code below commented out; it may be better to let the user decide
    % what to plot. Only make sure that no data are plotted.
    
    opt_tktd.preds = 1; % set to 1 only plot predictions from X0mat without data
    if ~(isempty(locL) && isempty(locR)) % if we have sub-lethal endpoints ...
        if opt_tktd.min == 1 % if user wants to see a dotted line for the control ...
            opt_tktd.addzero = 1; % set to 1 to always add a concentration zero to X0mat
            % Otherwise, the lowest plotted MF will be treated as the control
        end
    else % for survival, there is not much use since there is no control mortality
        opt_tktd.min     = 0; % set to 1 to show a dotted line for the control (lowest) treatment
    end
    
    make_scen(-5,-1); % remove all spline info for exposure profiles as we'll make a range of new ones
    
    Cw_C = Cw_tmp(2:end,2); % concentration vector of the exposure profile
    Cw_t = Cw_tmp(2:end,1); % time vector of the exposure profile
    scen = 0;
    % count and collect all MFs that need to be plotted
    MF_coll = [];
    for i_X = 1:length(ind_traits) % run through traits
        for i_F = 1:length(Feff) % run through effect levels
            if ~isnan(EPx{i_F}(i_X))
                scen = scen + 1; % count another plottable scenario
                MF_coll = cat(1,MF_coll,[EPx{i_F}(i_X) ind_traits(i_X) Feff(i_F)]); % collect this MF
            end
        end
    end
    MF_coll = sortrows(MF_coll,1); % sort based on MF
    
    X0mat   = [];
    for i = 1:scen % run through scenarios that need plotting
        
        switch MF_coll(i,2)
            case locS
                a = sprintf('%0.0f%% surv',100*MF_coll(i,3));
            case locL
                a = sprintf('%0.0f%% length',100*MF_coll(i,3));
            case locR
                a = sprintf('%0.0f%% repro',100*MF_coll(i,3));
            case loc_h
                a = sprintf('%0.0f%% hlthy',100*MF_coll(i,3));
        end
        
        make_scen(4,[1 i;Cw_t Cw_C * MF_coll(i,1)]); % create the globals to define the forcing function
        X0mat = cat(2,X0mat,[i;X0mat_tmp(2:end,1)]); % add a scenario to X0mat
        Label{i} = [a,' MF = ',num2str(round(MF_coll(i,1),3,'significant'))]; % create a label for plotting titles with the MF
    end
    Scenario = (1:scen)';
    Label    = Label';
    glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels
    
    glo.t = linspace(0,Cw_t(end),100);
    
    if ~isempty(opt_conf)
        opt_conf.set_zero = cat(2,backhaz,setzero); % parameter name(s) to set to zero (usually the background hazard hb)
        % This forces hb=0 (both in the best-fit parameter set and in the sample),
        % which is a good idea when simulating a long exposure profile (but not
        % when plotting a fit). However, that only works when we opt_conf is
        % not empty!
        if par_read == 1 % then we don't want to make CIs
            opt_conf = [];
        end
    end
    
    plot_tktd(par_plot,opt_tktd,opt_conf);
    % Note that in par_plot, the background hazard is set to zero
    
end

%% Return the globals to their original states

X0mat = X0mat_rem; 
glo   = glo_rem;   


function crit = calc_epx_sub(MF,t,par_out,Feff,X0_mat,Xctrl,locX,glo)

% Usage: crit = calc_epx_sub(MF,t,par_out,Feff,X0_mat,Xctrl,locX)
%
% This function calculates the criterion to be used by <fzero> to calculate
% the EPx. This function is generally useful for TKTD models, including
% GUTS and DEBtox. 
%
% Inputs:
% <MF>      the multiplication factor to try to see if it is the EPx
% <t>       the time vector for calculating the response of the traits
% <par_out> the optimised parameter set
% <Feff>    the fraction effect (x/100 in EPx)
% <X0_mat>  matrix with scenarios and initial values for the state variables
% <Xctrl>   value for the control trait (as EPx is relative to control)
% <locX>    location of trait of interest in Xout

glo.MF = MF; % modify global
Xtst2  = call_deri(t,par_out,X0_mat,glo); % response at conc. c (the identifier in X0mat)
crit   = (Xtst2(end,locX)/Xctrl)-(1-Feff); % zero when end value is x*100% effect


