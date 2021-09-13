function [LCx,LClo,LChi] = calc_lcx_lim_guts_red(par_plot,Tend,opt_lcx_lim,opt_conf)

% Usage: [LCx,LClo,LChi] = calc_lcx_lim_guts_red(par,Tend,opt_lcx_lim,opt_conf)
%
% Function to for fast calculation of LCx,t values, with confidence
% intervals. This function does the same as <calc_ecx>, but then faster,
% but ONLY works for standard GUTS models in this directory. That is to
% say: this function uses the complete analytical solutions for the LCx,t
% (IT) or the analytical solution for survival (SD), rather than calling
% call_deri. If you modify the model in <call_deri> or <simplefun>, this
% function could produce erroneous results (<calc_ecx> can still be used).
% 
% This function requires the model parameters to have their standard names.
% Since it depends on the specific names of parameters, this function is
% not placed in the engine but in the specific directory of the GUTS
% package that it applies to. Also, it assumes that initial damage (<Dw0>)
% is zero.
%
% <par_plot> parameter structure used for the best-fit line
% <Tend>     time point(s) at which to calculate LCx
% 
% Author     : Tjalling Jager 
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2

WRAP.glo  = glo;
WRAP.glo2 = glo2;

filenm    = glo.basenm;

% remember these globals as we are going to mess with them!
glo_rem = glo; % (make_scen changes it)

% read options from structure
Feff      = opt_lcx_lim.Feff; % effect level (>0 en <1), x/100 in LCx
make_plot = opt_lcx_lim.plot; % set to 0 to NOT make a plot of LCx vs time

if isempty(opt_conf)
    type_conf   = 0; % then we don't need CIs
    use_par_out = 0; % by default, set to off
else
    type_conf   = opt_conf.type; % use values from slice sampler (1), likelihood region(2) to make intervals
    use_par_out = opt_conf.use_par_out; % set to 1 to use par_out, as entered in this function, rather than from saved set
end
if type_conf > 0 || isempty(par_plot) % also if par_plot is not provided
    [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    if numel(rnd) == 1 % that means that no sample was found
        type_conf = 0; % no need to produce an error, just do LC50s without CI
    end
    if isempty(par_plot) % if no par structure was entered in this function ...
        par_plot = par; % simply use the one from the sample file
    else
        par_comp(par,par_plot,'hb') % compare the two parameter structures
    end
end

if ~isfield(glo,'fastslow') % does that field exist? 
    glo.fastslow = 'o'; % if not, set it 'off' ...
end
warning('off','backtrace')
if isfield(par_plot,'Di0') % just guessing ...
    warning('Using a parameter for initial damage requires careful consideration (it is NOT used here)')
    disp(' ')
end
if ~isempty(glo.names_sep)
    warning('You are using separate parameters per data set (with glo.names_sep)')
    warning('Only the FIRST set will be used for LCx (e.g., only f and not f1, f2, etc.')
    disp(' ')
end
warning('on','backtrace')

if numel(Feff) ~= 1
    error('This function cannot accept Feff as a vector; use one effect level at the time.')
end

if ~isfield(par_plot,'tag_fitted') % apparently, parameters have not been fitted
    warning('off','backtrace')
    warning('You did not fit any parameters, so LCx is based on the values in the parameter structure par that you entered.')
    warning('Any CIs are made from the saved set in the MAT file.')
    disp(' '), warning('on','backtrace')
end

make_scen(-5,-1); % remove all previously-defined scenario info (if any) from glo
% otherwise, predefined exposure scenarios may interfer

disp(' ')
disp('Calculate LCx values')
drawnow % clear the plotting buffer

%% Calculate LCx
% The sub-function calc_lc50 is used with <fzero> to find the LCx,t. Note
% that the sub-function creates an events sequence itself, assuming
% constant exposure (which is in the definition of LCx,t). Background
% hazard is set to zero.

Tend = sort(Tend); % make sure that time increases
% par.hb(1) = 0; % make background hazard zero (not needed, hb is not used anymore)

if Tend(1) == 0 % zero won't work as there is no effect there
    Tend(1) = []; % just remove the first entry
end
LCx = zeros(length(Tend),1); % initialise matrix to collect LCx,t
Tend = Tend(:); % this makes sure that LCx is a column vector

if glo.sel == 2 % for IT, we can use a direct calculation
    LCx = calc_lc50_it(Tend,par_plot,Feff,glo.fastslow);
elseif glo.sel == 1 && glo.fastslow == 'f' % for SD fast kinetics ...
    LCx = calc_lc50_sdf(Tend,par_plot,Feff); % we also have an analytical solution
else % for SD and mixed models, we use fzero
    LCinit = 2*par_plot.mw(1); % initial guess for LCx,t ... all should be between mw and max conc in data set
    LCx(1) = fzero(@calc_lc50_fzero,LCinit,[],Tend(1),par_plot,Feff,glo.sel,glo.fastslow); % find the first LCx,t
    if length(Tend) > 1
        for it = 2:length(Tend) % run through remaining time points for the LC50
            if glo.fastslow == 's' % for slow kinetics, mw is not a related to LC50 anymore
                LCx(it) = fzero(@calc_lc50_fzero,LCx(it-1),[],Tend(it),par_plot,Feff,glo.sel,glo.fastslow); % find the LCx,t 
                % initial value is value from previous T
            else
                LCx(it) = fzero(@calc_lc50_fzero,[par_plot.mw(1) 1.01*LCx(it-1)],[],Tend(it),par_plot,Feff,glo.sel,glo.fastslow); % find the LCx,t
                % make use of the fact that LC50s always decrease in time and is higher than mw!
            end
        end
    end
end

if make_plot == 1 % if required, make a plot!
    figh = make_fig(1,1); % create figure of correct size
    h_ax = gca;
    set(h_ax,'LineWidth',1,'FontSize',12) % adapt axis formatting
    hold on
    
    h_S = plot(Tend,LCx,'ko-','LineWidth',1,'MarkerFaceColor','y'); % remember handle for later bringing curve to the front
    
    xlabel('time (days)','Fontsize',12)
    ylabel(['LC',num2str(Feff*100),' (',glo.leglab2,')'],'Fontsize',12)
    xlim([0 Tend(end)])
    ylim([0 1.1*max(LCx)])
    
    if glo.saveplt > 0 && type_conf < 1 % if we want to save the plot and no CIs
        savenm = ['lcx_plot_',filenm];
        save_plot(figh,savenm);
    end   
end
drawnow % empty plot buffer if there's something in it

% par = par_plot; % unless we can take par from the saved file, use par_plot
% if type_conf > 0
%     [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
%     if numel(rnd) == 1 % that means that no sample was found
%         type_conf = -1; % no need to produce an error, just do LC50s without CI
%     end
%     par_comp(par,par_plot,'hb')
% end

diary('results.out') % collect output in the diary "results.out"
if type_conf < 1
    
    disp(['Results table for LCx,t (',glo.leglab2,')'])
    disp('==================================================================')
    for it = 1:length(Tend) % run through times requested
        fprintf('LC%1.0f (%4g days): %#5.4g \n',Feff*100,Tend(it),LCx(it))
    end
    disp('==================================================================')
    disp('No confidence intervals calculated')
    disp(' ')
    
    glo = glo_rem; % put the global back as it was
    diary off
    return % then we only calculate the best estimate, so we can return
end

%% Calculate confidence intervals
% Use the saved sample from the parameters-space explorer, and derive CIs
% from the results.

% NOTE: I think this function is so fast already that there is no need for
% using a parfor for the CIs.

if use_par_out == 1
    par = par_plot; % then we'll use the input par, rather than the saved one
    % Note: par_plot must already be structured in the main script, such
    % that the fitted parameters match the ones in the saved set, etc.
    % Note: when par_plot is NOT entered, it will have been made equal
    % to par from the saved set.
end

n_sets     = size(rnd,1); % number of samples in rnd
LCcoll     = zeros(n_sets,length(Tend)); % initialise matrix to collect LCx values
pmat       = packunpack(1,par,0,WRAP); % transform structure *from saved set* into a regular matrix
ind_fit    = (pmat(:,2)==1); % indices to fitted parameters
ind_logfit = (pmat(:,5)==0 & pmat(:,2)==1); % indices to pars on log scale that are also fitted!

f = waitbar(0,'Calculating CIs on LCx. Please wait.','Name','calc_lcx_lim.m');

for k = 1:n_sets % run through all sets in rnd
    
    waitbar(k/n_sets,f) % update waiting bar
    
    pmat(ind_fit,1) = rnd(k,:); % replace values in pmat with the k-th set in the sample
    % put parameters that need to be on log scale back on normal scale
    if sum(ind_logfit)>0
        pmat(ind_logfit,1) = 10.^(pmat(ind_logfit,1));
    end
    % Note: pmat is on normal scale here, but the sample in rnd contains
    % the value on a log scale, if a parameter is fitted on log scale.
    
    % Note: pmat is on normal scale here, but the sample in rnd contains
    % the value on a log scale, if a parameter is fitted on log scale.
    par_k = packunpack(2,0,pmat,WRAP); % transform parameter matrix into a structure
    % par_k.hb(1) = 0; % make background hazard zero NO NEED: hb is not used
    
    if glo.sel == 2 % for IT, we can use a direct calculation
        LCcoll(k,:) = calc_lc50_it(Tend,par_k,Feff,glo.fastslow);
    elseif glo.sel == 1 && glo.fastslow == 'f' % for SD fast kinetics ...
        LCcoll(k,:) = calc_lc50_sdf(Tend,par_k,Feff); % we also have an analytical solution
    else % for SD and mixed models, we use fzero
        if length(Tend) > 1
            for it = 1:length(Tend) % run through remaining time points for the LC50
                LCcoll(k,it) = fzero(@calc_lc50_fzero,LCx(it),[],Tend(it),par_k,Feff,glo.sel,glo.fastslow); % find the LCx,t
                % make use of the fact that LC50s will be close to best one
            end
        end
    end
       
end
close(f) % close the waiting bar
disp(' ')

% take min and max from the results of the sample
% this is done in this way to make sure that we don't mess up when either
% Tend or Teff is a single value (squeeze on 3-d matrix might return a
% swapped matrix).
disp(['Results table for LCx,t (',glo.leglab2,'), with 95% CI'])
disp('==================================================================')

LClo   = zeros(length(Tend),1);
LChi   = zeros(length(Tend),1);
for it = 1:length(Tend) % run through times requested
    switch type_conf % depending on the type of sample, take percentiles or min-max
        case 1
            LClo(it) = prctile(LCcoll(:,it),2.5,1); % 2.5 percentile of LCx
            LChi(it) = prctile(LCcoll(:,it),97.5,1); % 97.5 percentile of LCx
        case 2
            LClo(it) = min(LCcoll(:,it)); % take minimum
            LChi(it) = max(LCcoll(:,it)); % take maximum
        case 3
            LClo(it) = min(LCcoll(:,it)); % take minimum
            LChi(it) = max(LCcoll(:,it)); % take maximum
    end
    fprintf('LC%1.0f (%4g days): %#5.4g (%#5.4g-%#5.4g) \n',Feff*100,Tend(it),LCx(it),LClo(it),LChi(it))
end
disp('==================================================================')
switch type_conf
    case 1
        disp('Intervals: 95% credible (sample of posterior)')
    case 2
        disp('Intervals: 95% confidence (likelihood region)')
    case 3
        disp('Intervals: 95% confidence (parspace explorer)')
    otherwise
        disp('No confidence intervals calculated')
end
disp(' ')

if make_plot == 1 % if required, make a plot!
    figure(figh) % make the figure current, just to make sure we are plotting in the right place
    
    % Little trick to fill the area between the two curves, to
    % obtain a coloured confidence interval as a band.
    t2  = [Tend;flipud(Tend)]; % make a new time vector that is old one, plus flipped one
    Xin = [LClo;flipud(LChi)]; % do the same for the plot line, hi and lo
    fill(t2,Xin,'g','LineStyle','none','FaceAlpha',1) % and fill this object
    
    plot(Tend,LClo,'k:','LineWidth',1)
    plot(Tend,LChi,'k:','LineWidth',1)
    
    uistack(h_S,'top') % put model curve on top again, at the end of all plotting
    
    ylim([0 1.1*max(LChi)])
    switch type_conf
        case 1
            title('Intervals: 95% credible (sample of posterior','FontSize',10)
        case 2
            title('Intervals: 95% confidence (likelihood region)','FontSize',10)
    end
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['lcx_plot_',filenm];
        save_plot(figh,savenm);
    end
    
    drawnow
end

% put the globals back as they were
glo = glo_rem;

disp(['Time required: ' secs2hms(toc)])
diary off

%% ============== local function =========================================
function crit = calc_lc50_fzero(c,Tend,par,Feff,sel,fastslow)

% Usage: crit = calc_lc50_fzero(c,Tend,pmat,Feff,sel,fastslow)
%
% This function calculates the criterion to be used by <fzero> to calculate
% the LC50 for SD and mixed models. As LC50 is defines as 50% effect under
% constant exposure, we can use the analytical solution for survival here.
% <Tend> can be a vector.
%
% Inputs:
% <c>        the concentration to try to see if it is the LCx
% <Tend>     the time at which to calculate the LCx
% <par>      the optimised parameter set
% <Feff>     the fraction effect (x/100 in LCx)
% <sel>      SD or IT
% <fastslow> flag for fast/slow kinetics

kd   = par.kd(1); % dominant rate constant
mw   = par.mw(1); % median of threshold distribution
bw   = par.bw(1); % killing rate
Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution
beta = log(39)/log(Fs);       % shape parameter for logistic from Fs

switch sel
    case 1 % SD, this is the same analytical calculation as in DEBtool fomort
        S = 1;
        if fastslow == 's' % for slow kinetics, need to modify it as kd=0
            t0 = mw/c; % no-effect-time (note that mw is now a compound par!)
            f  = bw*mw*max(0,Tend-t0) - 0.5 * bw*c*max(0,Tend^2-t0^2);
            S  = min(1,exp(f)); % turn into survival probability
        elseif c > mw % if concentration to test is above the threshold ...
            t0 = -log(1-mw/c)/kd; % no-effect-time
            f  = (bw/kd)*max(0,exp(-kd*t0) - exp(-kd*Tend))*c - bw*(max(0,c-mw))*max(0,Tend - t0);
            S  = min(1,exp(f)); % turn into survival probability
        end
        
    case 3 % mixed model, this is the solution above plus some looping
        
        n          = 200; % number of slices from the threshold distribution
        Fs2        = 999^(1/beta); % fraction spread for 99.9% of the distribution
        z_range    = linspace(mw/(1.5*Fs2),mw*Fs2,n); % range of NECs to cover 99.9%
        prob_range = ((beta/mw)*(z_range/mw).^(beta-1)) ./ ((1+(z_range/mw).^beta).^2); % pdf for the log-logistic (Wikipedia)
        prob_range = prob_range / sum(prob_range); % normalise the densities to exactly one
        S          = 0; % initialise the survival probability over time with zeros
        for i = 1:n % run through the different thresholds
            surv = 1;
            if fastslow == 's' % for slow kinetics, need to modify it as kd=0
                t0   = mw/c; % no-effect-time (note that mw is now a compound par!)
                f    = bw*mw*max(0,Tend-t0) - 0.5 * bw*c*max(0,Tend^2-t0^2);
                surv = min(1,exp(f)); % turn into survival probability
            elseif c > mw % if concentration to test is above the threshold ...
                if fastslow == 'f' % for fast kinetics
                    surv = exp(-bw*(Cw-mw)*Tend); % shortcut
                else
                    t0   = -log(1-mw/c)/kd; % no-effect-time
                    f    = (bw/kd)*max(0,exp(-kd*t0) - exp(-kd*Tend))*c - bw*(max(0,c-mw))*max(0,Tend - t0);
                    surv = min(1,exp(f)); % turn into survival probability
                end
            end
            S     = S + surv * prob_range(i);  % add to S1, weighted for this NECs prob density
        end
        
end

crit = S-(1-Feff); % zero when end value is x*100% effect (for fzero to work on)

%% ============== local function =========================================

function LCx = calc_lc50_sdf(Tend,par,Feff)

% Usage: LCx = calc_lc50_sdf(Tend,par,Feff)
%
% This function calculates the LCx,t for SD when there is fast kinetics. As
% LC50 is defined as 50% effect under constant exposure, we can use an
% analytical solution here.  <Tend> can be a vector.

% extract parameters for the shortcut calculations
mw   = par.mw(1); % median of threshold distribution
bw   = par.bw(1); % killing rate

LCx = log(1/(1-Feff)) ./ (bw*Tend) + mw;

%% ============== local function =========================================

function LCx = calc_lc50_it(Tend,par,Feff,fastslow)

% Usage: LCx = calc_lc50_it(Tend,par,Feff,fastslow)
%
% This function calculates the LCx,t for IT. As LC50 is defined as 50%
% effect under constant exposure, we can use the analytical solution here.
% This also works when <Tend> is a vector.

% extract parameters for the shortcut calculations
mw   = par.mw(1); % median of threshold distribution
kd   = par.kd(1); % dominant rate constant
Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution
beta = log(39)/log(Fs);       % shape parameter for logistic from Fs

if fastslow == 's' % for slow kinetics, need to modify it as kd=0 (mw is then a compound par!)
    LCx = (mw./Tend) .* (Feff/(1-Feff))^(1/beta);
else
    LCx = (mw./(1-exp(-kd*Tend))) .* (Feff/(1-Feff))^(1/beta);
end

% %% ============== local function =========================================
% function crit = calc_lc50(c,Tend,par,Feff)
% 
% % Usage: crit = calc_lc50(c,Tend,pmat,Feff)
% %
% % This function calculates the criterion to be used by <fzero> to calculate
% % the LC50 for SD and the mixed model. This one is not used anymore, but
% % kept here for possible future use.
% 
% global glo % better not to make globals in sub-functions
% 
% Xtst = call_deri([0 Tend/2 Tend],par,[c X0'],glo); % survival at conc. c
% crit = Xtst(end,glo.locS)-(1-Feff); % zero when end value is x*100% effect

