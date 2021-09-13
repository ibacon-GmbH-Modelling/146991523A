function [LPx_out,LPlo,LPhi] = calc_lpx_lim_guts_full(par_plot,Tend,profname,opt_lcx_lim,opt_conf)

% Usage: [LPx_out,LPlo,LPhi] = calc_lpx_lim_guts_full(par,Tend,profname,opt_lcx_lim,opt_conf)
%
% Function for FAST calculation of LPx values, with confidence intervals.
% This function does the same as <calc_epx>, but then much faster, but ONLY
% works for standard GUTS models in this directory. For example: it will
% produce erroneous results for TK/damage models that use higher-order
% kinetics (e.g., saturation). The reason is that it makes use of the fact
% that, under first-order kinetics, the scaled damage level is proportional
% to the exposure concentration. Multiplying an entire exposure profile is
% thus equivalent to the multiplication of the damage levels over time by
% the same factor. We therefore only need to calculate <Dw> once for every
% parameter set. Furthermore, the survival probability is calculated here
% in this function. Therefore, it only works with the standard survival
% calculations for SD, IT and mixed. Further, furthermore, it requires that
% background hazard is zero when a parameter named <hb> is zero, and it
% requires the other parameters to have their standard names as well. Since
% it depends on the specific names of parameters, this function is not
% placed in the engine but in the specific directory of the GUTS package
% that it applies to. And finally, it requires initial damage (<Dw0>) to be
% zero. The standard function <calc_epx> is more generally applicable.
% 
% THIS IS THE VERSION BELONGING TO THE FULL GUTS MODEL
%
% LPx is the 'exposure-multiplication factor': the factor by which the
% exposure profile needs to be multiplied to obtain a certain percentage
% mortality after a specified duration. This only works when working with
% time-varying exposure profiles! (Otherwise, you might just as well
% calculate LCx,t and use that). 
% 
% <par_plot> parameter structure used for the best-fit line
% <Tend>     time point at which to calculate LPx (if left empty the end of the
%            exposure profile)
% <profname> either a file name with an exposure profile (time column and
%            concentration column) or the number of a scenario in the main script.
% 
% This function expects one effect level (in <opt_lcx_lim.Feff>) and one
% time point (<Tend>).
%
% Author     : Tjalling Jager
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 X0mat

WRAP.glo  = glo;
WRAP.glo2 = glo2;
% Note that glo will be changed in this function. That does not affect the
% functioning of WRAP, since WRAP is only used for packunpack here.

filenm    = glo.basenm;
batch_epx = 0; % don't use batch option for now

% rememeber these globals as we are going to mess with them!
glo_rem   = glo;
X0mat_rem = X0mat;

if ~isfield(glo,'fastslow') % does that field exist? 
    glo.fastslow = 'o'; % if not, set it 'off' ...
end
if numel(Tend)>1
    error('This function does not yet work for multiple time points. Call it several times with a new Tend.')
end

% read options from structure
Feff      = opt_lcx_lim.Feff; % effect level (>0 en <1), x/100 in LCx
make_plot = opt_lcx_lim.plot; % set to 0 to NOT make a plot of survival vs time at LPx
scen_plot = opt_lcx_lim.scen_plot; % set to 0 to NOT make a plot of the exposure profile
scen_type = opt_lcx_lim.scen_type; % type of definition for the exposure scenario (as used in make_scen), if called with a file name
notitle   = opt_lcx_lim.notitle; % set to 1 to suppress titles above ECx plots

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

if ~ismember(scen_type,[1 4])
    error('Only scenario type 1 or 4 are allowed at this point for LPx calculation, when calling this function with a named file with exposure scenario.')
end
warning('off','backtrace')
if ~isfield(par_plot,'tag_fitted') % apparently, parameters have not been fitted
    warning('You did not fit any parameters, so LPx is based on the values in the parameter structure par that you entered.')
    warning('Any CIs are made from the saved set in the MAT file.')
    disp(' ')
end
if ~isempty(glo.names_sep)
    warning('You are using separate parameters per data set (with glo.names_sep)')
    warning('Only the FIRST set will be used for LPx (e.g., only f and not f1, f2, etc.')
    disp(' ')
end
warning('on','backtrace')

if ischar(profname) % then we load a file with an exposure profile
    c = 1; % call the profile scenario 1
    glo.scen_plot = scen_plot; % this tells make_scen to make a plot or not
    make_scen(-5,-1); % remove all previously-defined scenario info (if any) from glo
    Cw_prof = load(profname,'-ascii'); % load the ringtest monitoring profile
    
    % Code below checks whether a pulsed scenario is used (with time points
    % occurring twice, which is used in openGUTS).
    if length(unique(Cw_prof(:,1))) ~= length(Cw_prof(:,1))
        % if we have a renewals scenario, with time points occurring
        % twice, extract a tiny bit from the first time point
        Cw_prof(diff(Cw_prof(:,1))==0,1) = Cw_prof(diff(Cw_prof(:,1))==0,1) - 0.01;
    end
    
    Cw_prof = [NaN 1; Cw_prof]; % add a first row with the scenario identifier 1
    make_scen(scen_type,Cw_prof); % prepare exposure scenario as linear forcing series
else % otherwise, assume the entered profname is a scenario number
    c = profname;
end
Tev = read_scen(-2,c,-1,glo); % use read_scen to derive actual exposure concentration
% the -2 lets read_scen know we need events, the -1 that this is also needed for splines!
if isempty(Tend)
    Tend = Tev(end,1); % by default take end of profile
elseif Tend > Tev(end,1)
    warning('off','backtrace')
    warning('You request a time point for LPx that is outside the specified exposure profile. Check the extrapolation!')
    disp(' '), warning('on','backtrace')
end

disp(' ')
disp(['Calculate exposure multiplication factors'])
drawnow % clear the plotting buffer

X0 = X0mat(2:end,1); % Note: by default, take the first column of X0mat for the 
% starting values of all state variable for the calculation of LPx.
% However, for the limited LPx calculation, we need some checks.
warning('off','backtrace')
if X0(glo.locD) ~= 0
    warning('Initial damage is forced to zero for this simplified LPx calculation (it is not zero in X0mat)!')
    disp(' ')
    X0(1+glo.locD) = 0;
end
if X0(glo.locS) ~= 1
    warning('Initial survival probablity is forced to one for this simplified LPx calculation (it is not one in X0mat)!')
    disp(' ')
    X0(1+glo.locS) = 1;
end
warning('on','backtrace')
if isfield(par_plot,'Ci0') % just guessing ...
    error('Using a parameter for initial concentration will lead to erroneous results)')
end
if isfield(par_plot,'Di0') % just guessing ...
    error('Using a parameter for initial damage will lead to erroneous results)')
end

%% Calculate LPx for the best parameter set

% First, construct a nice detailed time vector for the survival
% calculations (largest of 500 and 3 times the relevant part of the
% exposure profile). We need a lot of detail here, since we are using
% shortcuts for the LPx calculations that need a detailed vector for Dw.
min_t = 500;
t = Tev(Tev(:,1)<Tend,1); % only keep time points less than Tend
min_t = max(min_t,length(Tev(Tev(:,1)<t(end),1))*3); 
if length(t) < min_t % make sure there are at least min_t points
    t = unique([t;(linspace(0,Tend,min_t))']); 
    % this combines the points in the exposure scenario with a regular
    % vector along the exposure profile, and add Tend in there as well
end

% start by calculating Dw for multiplication factor 1
par_plot.hb(1) = 0; % make background hazard zero (needed for SD, to start initial value finder)
Xout      = call_deri(t,par_plot,[c;X0],glo); % survival and damage at scenario c
Di        = Xout(:,glo.locD); % scaled damage levels

if glo.sel == 2 % for IT, we can use a direct calculation
    
    Dim = max(Di); % find the maximum value for Dw over time
    LPx = calc_lpx_it(par_plot,Dim,Feff);
    
else % for SD, some more work is needed ... need to use fzero
    
    crit  = [1 Xout(end,glo.locS)-(1-Feff)]; % zero when end value is x*100% effect
    % positive when survival is higher than it should be, negative when
    % it is lower. For positive crit, we thus need to increase the
    % multiplication factor above 1. We start with the survival observed
    % with the normal result, so LPtry = 1.
    
    % we now step through multiplication factors, every time changing it by
    % a factor of 10, until we have an interval were the criterion changes
    % sign for each effect level (first going to higher factors, than to
    % lower, if needed)
    Di_mat = [t Di];
    LPtry = 1;
    
    while ~all(crit(end,2:end)<0) % continue until last row crit is negative
        LPtry = 10 * LPtry;
        crit1 = calc_lpx_sd(LPtry,par_plot,Feff,Di_mat,glo.sel);
        crit  = [crit ; LPtry crit1]; % add a new criterion at the bottom
    end
    LPtry = 1;
    while ~all(crit(1,2:end)>0) % continue until first row crit is positive
        LPtry = 0.10 * LPtry;
        crit1 = calc_lpx_sd(LPtry,par_plot,Feff,Di_mat,glo.sel);
        crit  = [LPtry crit1 ; crit]; % add a new criterion at the top
    end
    % crit becomes a matrix with multiplication factors and criteria
    % that will be used to provide an interval for fzero to search upon
    % (perhaps not really needed as fzero will find its way anyway;
    % seems to give some 5% time gain, so may be worth it).
    
    % calculate the LPx values for the best parameter set, for each effect level requested
    ind_init = [find(crit(:,2)>0,1,'last') find(crit(:,2)<0,1,'first')];
    LPx      = fzero(@calc_lpx_sd,crit(ind_init,1),[],par_plot,Feff,Di_mat,glo.sel); % find the LPx with fzero
    
end

% make a printout, before we get stuck into very long calculation of CIs ...
diary('results.out') % collect output in the diary "results.out"
disp('==================================================================')
fprintf('LP%1.0f (%4g days): %#5.4g \n',Feff*100,Tend,LPx)
disp('==================================================================')

LPx_out = LPx; % make the output variable, with the LPx value
LPlo    = -1; % these need to be defined if we stop prematurely
LPhi    = -1;

% par = par_plot; % unless we can take par from the saved file, use par_plot
% if type_conf > 0
%     [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
%     if numel(rnd) == 1 || isempty(rnd) % that means that no sample was found
%         type_conf = -1; % no need to produce an error, just do LC50s without CI
%     end
%     par_comp(par,par_plot,'hb') % compare the parameter structures in input and MAT file
%     par.hb(1) = 0; % make background hazard zero (needed for final plot)
% end

if make_plot == 1 % make a plot here already, without CIs
   
   glo.MF = LPx; % use the global to make read_scen know that we have an MF to apply
   % NOTE: the *input* par is used for the best curve!
   Xout = call_deri(t,par_plot,[c;X0],glo); % survival and damage at scenario c
   Di   = Xout(:,glo.locD); % scaled damage levels
   S    = Xout(:,glo.locS); % survival
   c_v  = read_scen(-3,c,t,glo); % use read_scen to derive exposure concentration vector
   
   figh = make_fig(1,2); % create a figure window of correct size
   hold on
 
   % subplot number 1 is for the exposure profile
   h_pl1 = subplot(1,2,1);
   hold on
   set(gca,'LineWidth',1,'FontSize',12) % adapt axis formatting
   xlabel(glo.xlab,'FontSize',12)
   if notitle == 0
       title(['Scenario: ',profname],'Fontsize',10,'Interpreter','none')
   end
      
   % This needs a little thought ... also how to deal with the option for
   % 'fast repair'. For now, don't do anything special when plotting.
   
%    if glo.fastslow == 's'
%        % for slow kinetics, damage is actually damage-time: the integrated
%        % damage over time. This makes it difficult to plot into the same
%        % graph as the exposure profile. Therefore, I use two y-axes.
%        yyaxis left
%        area(t,c_v,'LineStyle','none','FaceColor','b','FaceAlpha',0.25) % plot the exposure scenario as an area
%        ylabel(['concentration (',glo.leglab2,')'],'Fontsize',12)
%        xlim([0 Tend])
%        
%        yyaxis right
%        plot(t,Di,'k-','LineWidth',1)
%        ylabel(['damage-time (',glo.leglab2,' d)'],'Fontsize',12)
%        ylim([0 max(Di)*1.5])
%        xlim([0 Tend])
%        
%        legend({['profile x LP',num2str(100*Feff)],'damage-time'},'AutoUpdate','off')
%    else
       area(t,c_v,'LineStyle','none','FaceColor','b','FaceAlpha',0.25) % plot the exposure scenario as an area
       plot(t,Di,'k-','LineWidth',1)
       legend({['profile x LP',num2str(100*Feff)],'scaled damage'},'AutoUpdate','off')
       ylabel(['concentration (',glo.leglab2,')'],'Fontsize',12)
       xlim([0 Tend])
%    end
   
   % subplot number 2 is for the survival probability
   h_pl2 = subplot(1,2,2);
   hold on
   set(gca,'LineWidth',1,'FontSize',12) % adapt axis formatting
   
   h_S = plot(t,S,'k-','LineWidth',1); % remember handle for later bringing curve to the front
   
   xlabel(glo.xlab,'FontSize',12)
   ylabel('survival probability','Fontsize',12)
   ylim([0 1])
   xlim([0 Tend])
   title(['LP',num2str(Feff*100),': ',num2str(round(LPx,4,'significant')),' (',num2str(Tend),' d)'],'Fontsize',10)

   drawnow
   glo.MF = 1; % set it back to 1
   
   if glo.saveplt > 0 % if we want to save the plot, do it (will be overwritten by the one with CIs if needed)
       savenm = ['lpx_plot_',filenm];
       save_plot(figh,savenm);
   end
   
end

if type_conf < 1 % then no need to calculate CIs ... however we also skip the plot
    glo = glo_rem; % put the global back as it was
    diary off
    return % then we only calculate the best estimate, so we can return
end

%% Calculate confidence intervals
% Use the saved sample from the parameters-space explorer, and derive CIs
% from the results.

% Start/check parallel pool
if glo2.n_cores > 0
    poolobj = gcp('nocreate'); % get info on current pool, but don't create one just yet
    if isempty(poolobj) % if there is no parallel pool ...
        parpool('local',glo2.n_cores) % create a local one with specified number of cores
    end
end

if use_par_out == 1
    par = par_plot; % then we'll use the input par, rather than the saved one
    % Note: par_plot must already be structured in the main script, such
    % that the fitted parameters match the ones in the saved set, etc.
    % Note: when par_plot is NOT entered, it will have been made equal
    % to par from the saved set.
end
par.hb(1) = 0; % make background hazard zero in par (needed for final plot)

n_sets     = size(rnd,1); % number of sets in rnd
LPcoll     = zeros(n_sets,1); % initialise vector to collect LCx values for each set and effect level
pmat       = packunpack(1,par,0,WRAP); % transform structure *from saved set* into a regular matrix
ind_fit    = (pmat(:,2)==1); % indices to fitted parameters
ind_logfit = (pmat(:,5)==0 & pmat(:,2)==1); % indices to pars on log scale that are also fitted!
    
% if batch_epx == 0
%     f = waitbar(0,'Calculating CIs on LPx. Please wait.','Name','calc_lpx_lim.m');
% end

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
%     pmat(loc_zero,1) = 0; % make selected parameter(s) zero in each set of the sample!
%     % this has to be done after the tranformation to normal scale!
    pmat_coll{k} = pmat; % collect it in a huge cell array
end

if batch_epx == 0
    disp('Calculating CIs on LPx, please wait.')
end

% some trickery to get parfor running ...
glo_tmp    = glo;
sel        = glo.sel;
locD       = glo.locD;
glo_tmp.MF = 1; % set it back to 1 (was already done, but just to be sure ...)

parfor k = 1:n_sets % run through all sets in the sample
    
    par_k = packunpack(2,0,pmat_coll{k},WRAP); % transform parameter matrix into a structure
    % par_k.hb(1) = 0; % make background hazard zero NO NEED: hb is not used

    % calculate damage for this set from the sample
    Xout   = call_deri(t,par_k,[c;X0],glo_tmp); % survival and damage at scenario Tev
    Di     = Xout(:,locD); % scaled damage levels for this parameter set
    
    if sel == 2 % for IT, we can use a direct calculation
        
        Dim  = max(Di); % find the maximum value for Dw over time
        LPcoll(k) = calc_lpx_it(par_k,Dim,Feff);
                       
    else % for SD, some more work is needed ... need to use fzero
        
        % Again, make a crit matrix to help fzero with a range, starting at
        % the smallest value for LPx found for the best set (over the
        % various effect levels)
        Di_mat = [t Di]; % re-use the previously defined time vector
        
        LPtry = LPx; % start with factor for the basic run
        crit1 = calc_lpx_sd(LPtry,par_k,Feff,Di_mat,sel);
        crit  = [LPtry crit1]; % zero when end value is x*100% effect
        
        while ~all(crit(end,2:end)<0) % continue until last row all crit are negative
            LPtry = 10 * LPtry;
            crit1 = calc_lpx_sd(LPtry,par_k,Feff,Di_mat,sel);
            crit  = [crit ; LPtry crit1]; % add a new criterion at the bottom
        end
        LPtry = LPx; % restore LPtry to previous value
        while ~all(crit(1,2:end)>0) % continue until first row all crit are positive
            LPtry = 0.10 * LPtry;
            crit1 = calc_lpx_sd(LPtry,par_k,Feff,Di_mat,sel);
            crit = [LPtry crit1 ; crit]; % add a new criterion at the top
        end
        % crit becomes a matrix with multiplication factors and criteria
        % that will be used to provide an interval for fzero to search upon
        
        ind_init  = [find(crit(:,2)>0,1,'last') find(crit(:,2)<0,1,'first')];
        LPcoll(k) = fzero(@calc_lpx_sd,crit(ind_init,1),[],par_k,Feff,Di_mat,sel); % find the LPx with fzero
        
    end
    
end

% if batch_epx == 0
%     close(f) % close the waiting bar
% end
disp(' ')

% calculate CIs from the results of the sample
switch type_conf % depending on the type of sample, take percentiles or min-max
    case 1
        LPlo = prctile(LPcoll,2.5,1);  % 2.5 percentile of LPx
        LPhi = prctile(LPcoll,97.5,1); % 97.5 percentile of LPx
    case 2
        LPlo = min(LPcoll(:,1)); % take minimum
        LPhi = max(LPcoll(:,1)); % take maximum
    case 3
        LPlo = min(LPcoll(:,1)); % take minimum
        LPhi = max(LPcoll(:,1)); % take maximum
end

disp(['Results table for LPx (-) with 95% CI, scenario ',num2str(c)])
disp('==================================================================')
for ie = 1:length(Feff) % run through effect levels requested
    fprintf('LP%1.0f (%4g days): %#5.4g (%#5.4g-%#5.4g) \n',Feff(ie)*100,Tend,LPx(ie),LPlo(ie),LPhi(ie));
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
disp(['Time required: ' secs2hms(toc)])
diary off

if make_plot == 1
    
    glo.MF = LPx; % use the global to make read_scen know that we have an MF to apply
    
    opt_conf.samerr   = 0; % do not include sampling error in bounds for survival data (set samerr=1; requires statistics toolbox)
    if ~any(strcmp(opt_conf.set_zero,'hb')) % if hb is not specified to be set to zero yet 
        opt_conf.set_zero = cat(2,opt_conf.set_zero,'hb'); % parameter name to set to zero (usually the background hazard)
        % this adds hb to any existing parameter(s) that the user wants to set to zero
    end
    
    glo.t = t; % set time vector to exposure profile
    X0mat = [c;X0]; % and create a new X0mat for this profile
    out_conf = calc_conf(par,opt_conf); % calculate confidence intervals on model curves
    
    figure(figh); % make figure current again
    hold on
    
    % subplot number 1 is for the exposure profile
    subplot(h_pl1) % make this subplot current
    hold on
    
%     if glo.fastslow == 's'
%         % no need to plot CIs as there is no uncertainty on damage anymore
%     else
        % ylim([0 max(out_conf{2}{glo.locD})*1.5])
        plot(t,out_conf{1}{glo.locD},'k:','LineWidth',1)
        plot(t,out_conf{2}{glo.locD},'k:','LineWidth',1)
%     end
    
    % subplot number 2 is for the survival probability
    subplot(h_pl2) % make this subplot current
    hold on
    
    % Little trick to fill the area between the two curves, to
    % obtain a coloured confidence interval as a band.
    t2  = [t;flipud(t)]; % make a new time vector that is old one, plus flipped one
    Xin = [out_conf{1}{glo.locS};flipud(out_conf{2}{glo.locS})]; % do the same for the plot line, hi and lo
    fill(t2,Xin,'g','LineStyle','none','FaceAlpha',1) % and fill this object

    plot(t,out_conf{1}{glo.locS},'k:','LineWidth',1)
    plot(t,out_conf{2}{glo.locS},'k:','LineWidth',1)
   
    uistack(h_S,'top') % put model curve on top again, at the end of all plotting
    
    if glo.saveplt > 0 % if we want to save the plot and no CIs
        savenm = ['lpx_plot_',filenm];
        save_plot(figh,savenm);
    end
    
    glo.MF = 1; % set it back to 1
end

% put the globals back as they were
glo   = glo_rem;
X0mat = X0mat_rem;


%% ============== local function =========================================

function LPx = calc_lpx_it(par,Dwi,Feff)

mi   = par.mi(1);             % median of threshold distribution
Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution
beta = log(39)/log(Fs);       % shape parameter for logistic from Fs

LPx = (mi/Dwi) * (Feff./(1-Feff)).^(1/beta); % direct calculation of LPx

%% ============== local function =========================================

function crit = calc_lpx_sd(LPtry,par,Feff,Di_mat,sel)

% Usage: crit = calc_lpx_sd(LPtry,par,Feff,Di_mat,sel)
%
% This function calculates the criterion to be used by <fzero> to calculate
% the LPx. I make use of the fact that the scaled damage level is
% proportional to the exposure concentration. Multiplying an entire
% exposure profile is thus equivalent to the multiplication of the damage
% levels over time by the same factor. We therefore only need to calculate
% <Dw> once for every parameter set. This only works for standard
% first-order kinetics of damage (and internal concentration, in the full
% model)!
%
% Survival probability is calculated here; the code is largely copied from
% <call_deri>. This is used for SD and the GUTS-proper (IT and SD
% combined).
%
% Inputs:
% <LPtry>   the multiplication factor to see if it leads to x% effect
% <par>     the parameter set to use
% <Feff>    the fraction effect (x/100 in LPx)
% <Di_mat>  scaled damage over time, that we will multiply by <LPtry>
% <sel>     SD or IT(not used here) or GUTS-proper

mi   = par.mi(1);             % median of threshold distribution
bi   = par.bi(1);             % killing rate
Fs   = max(1+1e-6,par.Fs(1)); % fraction spread of the threshold distribution
beta = log(39)/log(Fs);       % shape parameter for logistic from Fs

% We want to calculate the multiplication factor to reach Feff at Tend. We
% already have the scaled damage levels Dw, which we can multiply with
% LPtry, and calculate the associated survival probability. This code is
% copied from call_deri.

t      = Di_mat(:,1);           % time vector in the data
Di     = LPtry * Di_mat(:,2);   % use multiplication on damage

switch sel
    case 1 % stochastic death
        haz    = bi * max(0,Di-mi);     % calculate hazard for each time point
        cumhaz = cumtrapz(t,haz);       % integrate the hazard rate numerically
        S      = min(1,exp(-1*cumhaz)); % calculate survival probability, incl. background
    
    case 3 % mixed model
        % This is a fast way for GUTS proper. Damage is calculated only
        % once, as it is the same for all individuals. Survival for
        % different thresholds is calculated below. 
        n          = 200; % number of slices from the threshold distribution
        Fs2        = 999^(1/beta); % fraction spread for 99.9% of the distribution
        z_range    = linspace(mi/(1.5*Fs2),mi*Fs2,n); % range of NECs to cover 99.9%
        prob_range = ((beta/mi)*(z_range/mi).^(beta-1)) ./ ((1+(z_range/mi).^beta).^2); % pdf for the log-logistic (Wikipedia)
        prob_range = prob_range / sum(prob_range); % normalise the densities to exactly one
        S          = zeros(length(t),1); % initialise the survival probability over time with zeros
        for i = 1:n % run through the different thresholds
            haz    = bi * max(0,Di-z_range(i));  % calculate hazard for each NEC
            cumhaz = cumtrapz(t,haz);            % integrate the hazard rate numerically
            surv   = min(1,exp(-1*cumhaz));      % calculate survival probability
            S      = S + surv * prob_range(i);  % add to S1, weighted for this NECs prob density
        end
end

crit   = S(end)-(1-Feff); % zero when end value is x*100% effect


