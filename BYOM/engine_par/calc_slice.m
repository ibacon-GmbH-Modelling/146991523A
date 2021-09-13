function [rnd,coll_rndsort] = calc_slice(par,nrs,opt_slice)

% Usage: [rnd,coll_rndsort] = calc_slice(par,nrs,opt_slice)
%
% Calculation of posterior probability (optional informative priors in
% script file) using MCMC (the <slicesample> algorithm in the statistics
% toolbox). Marginal densities are calculated from the sample using the
% kernel density estimation in the statistics toolbox (ksdensity).
% Confidence intervals are calculated from the densities by finding the
% parameter values with the highest densities that together yield 95% of
% the total area under the curve. This means that 'islands' of likelihood
% are also captured in the confidence intervals. Additionally, quantiles
% are provided.
%
% <nrs> is the number of samples to be drawn. This function saves the
% random sample into a file with as name: filenm_MC.MAT (<filenm> is the
% name of the mydata file that called the slice sampler). This saved file
% can be used to calculate prediction intervals. Optional output is the
% random sample and the CIs of the parameters (as 0.025-0.975 quantiles).
%
% Options can be set in a structure as third argument: <opt_slice>. This is
% predefined with defaults in <prelim_checks>.
%
% Author     : Tjalling Jager
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2
global DATA W DATAx Wx X0mat

WRAP.DATA  = DATA;
WRAP.W     = W;
WRAP.DATAx = DATAx;
WRAP.Wx    = Wx;
WRAP.X0mat = X0mat;
WRAP.glo   = glo;
WRAP.glo2  = glo2;

names = glo2.names;
pri   = glo2.pri;

% read options from structure
thin         = max(1,opt_slice.thin); % minimum is 1; thinning of the sample (keep one in every 'thin' samples)
burn         = opt_slice.burn;        % number of burn-in samples (0 is no burn in)
slwidth      = opt_slice.slwidth;     % initial width of the slice (Matlab default is 10) 
alllog       = opt_slice.alllog;      % put all parameters on log-scale before taking the sample
testing      = opt_slice.testing;     % make additional tests on the sample: moving average and autocorrelation
use_subplots = opt_slice.subplt;      % create single plot with subplots (1) or many individual plots (0)
filenm       = glo.basenm;

% By default, this file will use the name of the base script to load a MAT
% file. However, we also want to use saved MAT files for predictions in new
% script files. The glo.mat_nm can be used for that.
if isfield(glo,'mat_nm')
    filenm_load = glo.mat_nm;
else
    filenm_load = glo.basenm;
end

disp(' ')
warning('off','backtrace')
if exist('slicesample','file')~=2 % only when slicesample exists as an m-file in the path
    error('The function calc_slice needs the statistics toolbox installed ...')
    % I need slicesample and ksdensity, and copying them is not a good idea ...
end
if ~isfield(par,'tag_fitted') % apparently, parameters have not been fitted
    warning('Parameters have not been fitted, so results may not be very meaningful!')
    disp(' ')
end
warning('on','backtrace')

drawnow % plot the last plot in the plot buffer before starting the analysis

names_tmp = names; % work with a copy of names, as the saved set may have different names!
% if the saved set has different names, this likely will produce an error
% elsewhere anyway.

pmat = packunpack(1,par,0,WRAP); % transform par structure into a regular matrix

if alllog == 1 % if we want the slice sampler to operate on log parameters
    % Note: this does not affect the priors: they are defined on normal
    % scale, and continue to operate on normal scale. Their plotting is,
    % however, affected, as log-parameters are plotted on log-scale.
    pmat(:,5) = 0; % put all parameters on log acale
    pmat(:,3) = max(pmat(:,3),1e-20); % make sure that the lower bound is not zero
    pmat(:,4) = min(pmat(:,4),1e+20); % make sure that the upper bound is not infinite
    pmat(:,1) = max(pmat(:,1),1e-20); % make sure that the start value is not zero
end

pmat_tmp = pmat; % use this to compare to the saved one, if re-using a saved set

%% If re-using a saved data set, load it first

if nrs == -1
    if exist([filenm_load,'_MC.mat'],'file') ~= 2
        error('There is no MCMC sample saved, so run calc_slice again with positive number of samples')
    end
    disp('Loading and running the previously saved set. Also the parameter matrix is taken from the saved set.')
    load([filenm_load,'_MC'],'rnd','par','par_sel') % load the random sample from the last MCMC run
    % this loads the random sample in rnd, and the parameter matrix par and the
    % selection matrix par_sel.
    
    % rnd is the full sample, with a log-likelihood value in the last column
%     acc = sortrows(rnd,size(rnd,2));
%     acc = acc(acc(:,end)>=critloglik,1:end-1); % these are the accepted parameters (within the 95% credible region)
    rnd = rnd(:,1:end-1); % remove likelihood column from rnd (not used now)
    
    % repeat some of the previous actions, as par is now replaced with the
    % matrix from the saved set!
    pmat = packunpack(1,par,0,WRAP);  % transform globals structure into a regular matrix
    % Note that log-parameters are on normal scale in par, and thus also in pmat here
    names_tmp  = fieldnames(par); % extract all field names of par (as saved)
    
    if ~isequal(pmat,pmat_tmp) % if the saved set differed from the par as entered ...
        disp(' '), warning('off','backtrace')
        if isequal(pmat(:,1:2),pmat_tmp(:,1:2)) % ah, it's just the log settings or boundaries
            warning('The log settings and/or boundaries of parameters in the saved set differs from that in the workspace at the moment. This may not hinder the analysis.')
        else
            warning('The saved parameter matrix appears to be (slightly) different from the one entered when calling calc_slice. The one from the saved set is used!')
        end
        warning('on','backtrace'), disp(' ')
        fprintf('Parameter values from saved set \n');
        fprintf('=================================================================================\n');
        nfields = length(names);
        for i = 1:nfields % display all parameters on screen
            if pmat(i,5) == 0
                fprintf('%-6s %10.4g (fit: %1.0f) bounds: %6.4g - %6.4g log-scale \n',names_tmp{i},pmat(i,1),pmat(i,2),pmat(i,3),pmat(i,4))
            else
                fprintf('%-6s %10.4g (fit: %1.0f) bounds: %6.4g - %6.4g \n',names_tmp{i},pmat(i,1),pmat(i,2),pmat(i,3),pmat(i,4))
            end
        end
        fprintf('=================================================================================\n');
        fprintf('  \n');
    end
end

% put parameters that need to be on log scale on log10 scale; this is what transfer wants
pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));
% do NOT modify the min max ranges, as transfer needs them on normal scale

p_index = find(pmat(:,2)==1);
parshat = pmat(p_index,1); % Only the parameters that were fitted (log-transformed if needed)
par_sel = p_index; % this one will be saved with the MAT file later

%% Run slicesample.m (statistics toolbox) for the MCMC sample

if nrs ~= -1
    if alllog == 1 % when me modified the parameters, change par (which is saved in the MAT file)
        par = packunpack(2,0,pmat_tmp,WRAP);  % transform regular matrix into globals structure
        % Note that log-parameters are on normal scale in par
    end
    disp('Calculating sample from posterior density ... please be patient.')
    
    opt_test = 0; % 0) default, 1) option to use modified slice sampler
    
    if opt_test == 0
        rnd = slicesample(parshat,nrs,'logpdf',@(pars) -1 * transfer(pars,pmat,WRAP),'thin',thin,'burnin',burn,'width',slwidth);
        % transfer gives the min log likelihood, here we need the log likelihood, so multiply by -1
    else % modifed slicesample, provided in BYOM, which can output the MLL and parameter sets
        % this option is, at the moment, restricted to use by Tjalling
        rnd = slicesample_byom(parshat,nrs,pmat,WRAP,thin,burn,slwidth,0);
    end        
end

if length(p_index) == 1 % when we fit only 1 parameter ...
    use_subplots = 0; % turn of the subplots (no use for them)
end
if testing == 1 % more detailed testing of the sample
    testslice(rnd,names_tmp(p_index))
end

%% Plot marginal distributions and calculate credible intervals

corr_rnd = corrcoef(rnd); % correlation matrix from the sample
rndsort  = sort(rnd); % sort the MCMC sample (sorts each column in ascending order)

diary (glo.diary) % collect output in the diary "results.out"
disp(' ')
disp('Median and credible intervals (95% highest posterior density, 2.5-97.5% quantiles)')
disp('======================================================================================')

% % plotting the cumulative posterior in raw form
% cdf_x = linspace(0,1,size(rnd,1)); % for cumulative density plots, if needed
% for parnr = 1:length(p_index), % run through all fitted parameters
%     figure
%     hold on
%     if pmat(p_index(parnr),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
%         xlabel(['parameter ',names_tmp{p_index(parnr)},' (log-scale)'],'FontSize',12,'FontName','Arial')
%     else
%         xlabel(['parameter ',names_tmp{p_index(parnr)}],'FontSize',12,'FontName','Arial')
%     end
%     ylabel('cumulative density','FontSize',12,'FontName','Arial')
%     set(gca,'LineWidth',1,'FontSize',12,'FontName','Arial')
%     plot(rndsort(:,parnr),cdf_x,'k-') % plot cumulative density, if needed
% end

if use_subplots == 1
    [figh,ft] = make_fig(length(p_index),length(p_index)); % create figure of correct size
    hold on
end

maxrnd = zeros(length(p_index),1); % pre-define with zeros
minrnd = zeros(length(p_index),1); % pre-define with zeros
coll_rndsort = zeros(length(p_index),3); % pre-define with zeros
coll_rndsort(:,3) = median(rndsort)'; % put the median in last column of coll_rndsort

for parnr = 1:length(p_index) % run through all fitted parameters
    
    if use_subplots == 0
        make_fig(1,1); % create figure of correct size
        set(gca,'LineWidth',1,ft.name,ft.ticks)
        title(['Called from: ',filenm, ' (',date,')'],'Interpreter','none',ft.name,ft.title)
        if pmat(p_index(parnr),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
            xlabel(['parameter ',names_tmp{p_index(parnr)},' (log-scale)'],ft.name,ft.label)
        else
            xlabel(['parameter ',names_tmp{p_index(parnr)}],ft.name,ft.label)
        end
        ylabel('estimated density',ft.name,ft.label)
    else
        g = subplot(length(p_index),length(p_index),(parnr-1)*length(p_index)+parnr);
        set(gca,'LineWidth',1,ft.name,ft.ticks)
        p = get(g,'position');
        p([3 4]) = p([3 4])*1.25; % add 10 percent to width and height
        if length(p_index) == 2 % when plotting just 2 pars, need to move them down a bit
            p(2) = p(2)-0.04;
        end
        set(g, 'position', p);
        ax = gca; % handle to current axis system

        if parnr == 1 % for first plot, give it a y-axis
            if pmat(p_index(parnr),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                ylabel([names_tmp{p_index(parnr)},' (log)'],ft.name,ft.label)
            else
                ylabel([names_tmp{p_index(parnr)}],ft.name,ft.label)
            end
            set(ax,'XTickLabel',[]); % this works with old versions as well
        elseif parnr == length(p_index) % for last plot, give it an x-axis
            if pmat(p_index(parnr),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                xlabel([names_tmp{p_index(parnr)},' (log)'],ft.name,ft.label)
            else
                xlabel([names_tmp{p_index(parnr)}],ft.name,ft.label)
            end
            set(ax,'YTickLabel',[]); % 
        else
            set(ax,'XTickLabel',[]); % 
            set(ax,'YTickLabel',[]); % 
        end
        
    end
    hold on
     
    % [ind1,ind2] = ind2sub(size(par_sel),p_index(parnr)); % index to the fitted parameter in 2 dims

    if pmat(p_index(parnr),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
        supp_par = log10(pmat(p_index(parnr),[3 4])); % support for kernel density is min-max on log scale
    else
        supp_par = pmat(p_index(parnr),[3 4]); % support for kernel density is min-max on normal scale
    end
    
    [f,x,~] = ksdensity(rndsort(:,parnr),'support',supp_par); % use kernel density estimation
    % Note: limit support for the posterior to the min-max range provided in the byom script
    plot(x,f,'k-'); % plot the estimated marginal probability density for the parameter
    
    if isfield(pri,names_tmp{p_index(parnr)}) % is there a prior for this one?
        % If no prior is provided for a parameter, a uniform one (on normal
        % scale) is taken. This is NOT plotted (and would be a bit
        % complicated for log-scale parameters).
        
        ptmp = []; % clear the temporary structure
        
        % find a reasonable range to plot the prior on
        a = eval(['pri.',names_tmp{p_index(parnr)}]); % extract prior info
        switch a(1)
            case 1
                pribnds = [0 1];
            case 2
                pribnds = [a(2) a(3)];
            case 3
                pribnds = [a(2)-3*a(3) a(2)+3*a(3)];
            case 4
                pribnds = exp([a(2)-3*a(3) a(2)+3*a(3)]);
        end
        if pmat(p_index(parnr),5) == 0 % if parameter is fitted on log-scale
            pribnds(1) = max(pribnds(1),10^supp_par(1));
            pribnds = log10(pribnds); % put them on log scale
        end
        
        % set plotting range to a reasonable range to catch the entire
        % posterior and a relevant part of the prior
        % rangep = linspace(min(x(1),max(pribnds(1),x(1)-abs(x(1)/2))),max(x(end),min(pribnds(2),x(end)+abs(x(end)*2))),500); % range of values for plotting prior
        
        rangep = linspace(pribnds(1),pribnds(2),500); % range of values for plotting prior
        % I now take the (almost) full range so that I can scale when the
        % parameter is on log scale. This should help to get a properly
        % scaled plot of the prior.
               
        priprob = zeros(length(rangep),1);
        for ip = 1:length(rangep)
            if pmat(p_index(parnr),5) == 0 % if parameter is fitted on log-scale
                eval(['ptmp.',names_tmp{p_index(parnr)},'=10.^rangep(ip);']);
            else
                eval(['ptmp.',names_tmp{p_index(parnr)},'=rangep(ip);']);
            end
            priprob(ip) = exp(calc_prior(ptmp,glo2)); % call function to calculate prior prob (exponential as priprob gives log-prob values).
            
        end
        
        if pmat(p_index(parnr),5) == 0 % if parameter is fitted on log-scale
            priprob = priprob / sum(priprob*(rangep(2)-rangep(1))); % rescale the prior for proper plotting
        end
                
        plot(rangep,priprob,'r-')
        if use_subplots == 0
            xlim([rangep(1) rangep(end)])
        end
         
    end
    
    % Next, find the critical density for the 95% highest-prob. region, using fzero
    stepsize = x(2)-x(1);
    clear coll_int
    critdens = fzero(@(testdens) sum(f(f>=testdens))*stepsize-0.95,max(f)/2); % this is the prob density that leads to 95% of the values
    plot([min(x) max(x)],[critdens critdens],'k:') % plot this cut off prob density as a horizontal line
    
    if isnan(critdens) % in extreme cases, a critical density cannot be found
        coll_int = [nan nan];
    else
        % Now calculate the 95% highest-prob. boundaries
        clear coll_int
        testrange = f-critdens;
        direction = -1; % f-critdens starts as negative (we always first find the start of an interval)
        nr_int = 0;
        for i = 1:length(x)
            if direction == -1 && testrange(i) > 0 % we found the start of an interval
                nr_int = nr_int+1;
                if i > 1
                    coll_int(nr_int,1) = mean([x(i) x(i-1)]); % take average as start of interval
                else
                    coll_int(nr_int,1) = x(1); % take first as start of interval
                end
                direction = 1; % change the direction to spot end of interval
            end
            if direction == 1 && testrange(i) < 0 % we found the end of an interval
                if i<length(x)
                    coll_int(nr_int,2) = mean([x(i) x(i-1)]); % take average as end of interval
                else
                    coll_int(nr_int,2) = x(end); % take last as end of interval
                end
                direction = -1; % change the direction to spot start of new interval
            end
            if direction == 1 && testrange(i) > 0 && i == length(x) % we found no end of an interval at the last x value
                coll_int(nr_int,2) = x(end); % take last as end of interval
            end
        end
    end
    
    % Calculate the quantiles

    %     coll_rndsort(parnr,:) = quantile(rndsort(:,parnr),[0.025 0.975]); % numbers for quantiles
    %     % plot broken lines to indicate quantiles
    %     plot([coll_rndsort(parnr,1) coll_rndsort(parnr,1)],[0 interp1(x,f,coll_rndsort(parnr,1))],'k--')
    %     plot([coll_rndsort(parnr,2) coll_rndsort(parnr,2)],[0 interp1(x,f,coll_rndsort(parnr,2))],'k--')
    
    % using the kernel densities for the quantiles is better for small samples
    if any(f<1e-20) % sometimes there are zeros or very low numbers in there!
        f = f+1e-10; % add a small amount should prevent error...
    end
    coll_rndsort(parnr,1) = interp1(cumtrapz(x,f),x,0.025);
    coll_rndsort(parnr,2) = interp1(cumtrapz(x,f),x,0.975);
    % If this still yields an error: try with more samples or use the commented-out code above
    
    % plot broken lines to indicate quantiles
    plot([coll_rndsort(parnr,1) coll_rndsort(parnr,1)],[0 interp1(x,f,coll_rndsort(parnr,1))],'k--')
    plot([coll_rndsort(parnr,2) coll_rndsort(parnr,2)],[0 interp1(x,f,coll_rndsort(parnr,2))],'k--')
    
    if use_subplots == 0
        axis tight % axis as close as possible
    else % remember min and max values so that sample plots can also be scaled in the same way
        maxrnd(parnr) = max(max(x),max(rndsort(:,parnr)));
        minrnd(parnr) = min(min(x),min(rndsort(:,parnr)));
        xlim([minrnd(parnr) maxrnd(parnr)])
    end
    
    if glo.saveplt > 0 && use_subplots == 0 % if we want to save the plot, and don't make subplots
        savenm = ['slice_',filenm,'_par_',names_tmp{p_index(parnr)}];%
        save_plot(gcf,savenm);
    end
    
    if pmat(p_index(parnr),5) == 0 % parameter is on log scale
        coll_rndsort(parnr,:) = 10.^coll_rndsort(parnr,:); % for screen, make output normal scale
        coll_int = 10.^coll_int;
    end
    
    if nr_int == 1 % no problems, it is a single interval
        % fprintf('%-5s (%1.0f) interval: %10.4g - %10.4g quant.: %10.4g - %10.4g\n',PARS{ind1,1},ind2,coll_int,coll_rndsort(parnr,:))
        fprintf('%-10s %10.4g hdr: %10.4g - %10.4g quant: %10.4g - %10.4g\n',names_tmp{p_index(parnr)},coll_rndsort(parnr,3),coll_int,coll_rndsort(parnr,1:2))
    else
        disp('The credible interval is a broken set (check density plots to check these figures)')
        fprintf('%-10s %10.4g                      quant: %10.4g - %10.4g\n',names_tmp{p_index(parnr)},coll_rndsort(parnr,3),coll_rndsort(parnr,1:2))
        for i_broken = 1:nr_int
            fprintf('           hdr: %10.4g - %10.4g \n',coll_int(i_broken,:))
        end
    end
    
end

%% Find the joint 95% highest-prob. region
% Here, I calculate the log-likelihood for each parameter set to be able to
% construct a joint, highest-probability, credible region! This is also
% saved to the file. NOT USED ANYMORE.

if nrs ~= -1 % don't calculate acc and save again when we use the saved sample anyway
    minloglik_rnd = zeros(size(rnd,1),1); % pre-define with zeros
    
%     for i = 1:nrs % run through all samples
%         minloglik_rnd(i) = transfer(rnd(i,:),pmat,WRAP); % use transfer to obtain minloglikelihood!
%     end % pmat includes the log parameters on log scale, as transfer needs
    
% NOTE: I removed the option to have a limited set with the highest
% posterior density. I think it makes little sense to calculate it, and the
% calculation of the minloglik may be time consuming. In the future, I may
% implement a different way to use the sample (in a frequentist manner) so
% I leave the option in, commented out.

    rnd  = [rnd minloglik_rnd]; % add the log-liklihood to the random sample
%     rnd  = sortrows(rnd,size(rnd,2)); % sort based on likelihood
%     acc  = rnd(1:floor(0.95*nrs),1:end-1); % take the 95% best ones
    % NOTE: don't sort rnd itself before saving, as that precludes the
    % possibility to use test_slice on a saved sample!
    
    % save the sample in a MAT-file with the name of the script, with
    % _MC at the end. Load can be used to retrieve it.
    % also, save par (total parameter matrix) and par_sel (selection
    % matrix) to make plots without redoing the fit.
    save([filenm,'_MC'],'rnd','par','par_sel')
    
    % If we're doing a DEBtox analysis, it is a good idea to save the
    % brood-pouch delay with the sample. 
    if isfield(glo,'Tbp')
        Tbp = glo.Tbp;
        save([glo.basenm,'_MC'],'Tbp','-append')
    end
    % If we have multiple versions of the same parameter, also save that info.
    if isfield(glo,'names_sep')
        names_sep = glo.names_sep;
        save([glo.basenm,'_MC'],'names_sep','-append')
    end

    rnd = rnd(:,1:end-1); % remove the likelihood column again
end

%% Plot the sample for each parameter combination

if use_subplots == 1

    for parnr1 = 1:length(p_index)-1 % go through columns
        for parnr2 = parnr1+1:length(p_index) % go through rows
            figure(figh) % make sure that the multiplot is the current plot
            g = subplot(length(p_index),length(p_index),(parnr2-1)*length(p_index)+parnr1);
            set(gca,'LineWidth',1,ft.name,ft.ticks)
            
            p = get(g,'position');
            p([3 4]) = p([3 4])*1.25; % add to width and height
            if length(p_index) == 2 % when plotting just 2 pars, need to move them down a bit
                p(2) = p(2)-0.04;
            end
            set(g, 'position', p);
            ax = gca;
            hold on
            
            mu     = mean([rnd(:,parnr1) rnd(:,parnr2)]);
            covmat = cov([rnd(:,parnr1) rnd(:,parnr2)]);
            h = error_ellipse(covmat,mu,'conf',0.95); % calculate and plot an error ellipse
            set(h,'LineWidth',1.5); % make line for error ellipse thicker
            
            plot(rnd(:,parnr1),rnd(:,parnr2),'ko','MarkerFaceColor','w')
            % plot(acc(:,parnr1),acc(:,parnr2),'ko','MarkerFaceColor','y')
            plot(parshat(parnr1),parshat(parnr2),'ko','MarkerFaceColor','r','MarkerSize',8)
            
            axis([minrnd(parnr1) maxrnd(parnr1) minrnd(parnr2) maxrnd(parnr2)])
            
            if parnr1 == 1 % we are in the first column
                if pmat(p_index(parnr2),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                    ylabel([names_tmp{p_index(parnr2)},' (log)'],ft.name,ft.label)
                else
                    ylabel([names_tmp{p_index(parnr2)}],ft.name,ft.label)
                end
                if parnr2<length(p_index) % if we are not on the last row, remove axis numbers
                    set(ax,'XTickLabel',[]); % this works with old versions as well
                else
                    if pmat(p_index(parnr1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                        xlabel([names_tmp{p_index(parnr1)},' (log)'],ft.name,ft.label)
                    else
                        xlabel([names_tmp{p_index(parnr1)}],ft.name,ft.label)
                    end
                end
            elseif parnr2 == length(p_index) % we are in the last row
                if pmat(p_index(parnr1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                    xlabel([names_tmp{p_index(parnr1)},' (log)'],ft.name,ft.label)
                else
                    xlabel([names_tmp{p_index(parnr1)}],ft.name,ft.label)
                end
                if parnr1>1 % if we are not on the first column, remove axis numbers
                    set(ax,'YTickLabel',[]); % this works with old versions as well
                end
            else
                set(ax,'YTickLabel',[]); % this works with old versions as well
                set(ax,'XTickLabel',[]); % this works with old versions as well
            end
            
        end
    end
    
    % create a legend in top-right sub-plot
    h_annot = subplot(length(p_index),length(p_index),length(p_index)); % make an extra subplot, and remember the handle
    set(gca,'LineWidth',1,ft.name,ft.ticks)
    p = get(h_annot,'position');
    p([3 4]) = p([3 4])*1.2; % add to width and height
    if length(p_index) == 2 % when plotting just 2 pars, need to move them down a bit
        p(2) = p(2)-0.04;
    end
    h1 = gca; % remember the current axis
    dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
    hold on
    
    % plot fake data points, with same symbols as in regular plots
    plot(0,0,'ko','MarkerFaceColor','w')
    plot(0,0,'ko','MarkerFaceColor','r','MarkerSize',8)
    
    plot([0 0],[0.1 0.1],'k-')
    plot([0 0],[0.1 0.1],'k--')
    plot([0 0],[0.1 0.1],'k:')
    
    % define legend entries
    L = {'sample from posterior','best-fitting set','kernel-density estimate','quartiles on parameters','cutoff max. posterior density'};
    h_leg_tot = legend(L); % create a legend with all entries
    set(h_leg_tot,ft.name,ft.legend); % use standard font formats
    
    xlim([1 10]) % make sure the 'data' are off screen, so invisible
    dimleg = h_leg_tot.Position; % position of legend made
    % move it to the new subplot and set it to top-left of panel
    h_leg_tot.Position = [dim(1) dim(2)+dim(4)-dimleg(4) dimleg(3:4)];
    h1.Position = [dim(1) dim(2)+dim(4)-dimleg(4) 0.01 0.01];
    h1.Visible = 'off'; % hide axes
    
    if ~isfield(glo,'saveplt_notitle') || glo.saveplt_notitle ~= 1
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Bayesian sample from posterior. Using file: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
    end
    if glo.saveplt > 0 % if we want to save the plot, and make subplots
        savenm = ['slice_multi_',filenm];%
        save_plot(gcf,savenm);
    end
    
%     % Include a text box with the results
%     h_annot = subplot(length(p_index),length(p_index),2); % make an extra subplot, and remember the handle
%     set(gca,'LineWidth',1,'FontSize',10,'FontName','Arial')
%     p = get(h_annot,'position');
%     p([3 4]) = p([3 4])*1.2; % add to width and height
%     if length(p_index) == 2 % when plotting just 2 pars, need to move them down a bit
%       p(2) = p(2)-0.04;
%     end
%     set(h_annot, 'position', p); % this minimises the white space between the plots
%     h1 = gca; % remember the current axis
%     dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
%     str_annot = []; % clear the string if it already exists, and fill it step by step
%     
%     str_annot{1} = 'Quantiles of posterior';
%     str_annot{2} = '=========================';
%     for parnr = 1:length(p_index) % run through all fitted parameters
%         str_annot{2+parnr} = sprintf('%-4s %1.3g - %1.3g',names_tmp{p_index(parnr)},coll_rndsort(parnr,1),coll_rndsort(parnr,2));
%     end
%     str_annot{end+1} = '=========================';
%     
%     % make a textbox and place it in the additional subplot
%     annotation('textbox',dim,'String',str_annot,'FitBoxToText','on','BackGroundColor','w','FontName','Courier','Interpreter','none');
%     delete(h1) % throw away the axes of the plot, the text box is enough
    
end

%% Finish printing on screen 

disp('======================================================================================')
disp(' ')
disp('Correlation matrix from the MCMC sample')
disp('======================================================================================')
fprintf('%5s ','  ')
for i = 1:length(p_index)
    % [inda,indb]=ind2sub(size(par_sel),p_index(i)); % finds the subscripts from the indices of this parameter
    fprintf('%11s ',names_tmp{p_index(i)})
end
fprintf('\n')
for i = 1:length(p_index)
    % [inda,indb]=ind2sub(size(par_sel),p_index(i)); % finds the subscripts from the indices of this parameter
    fprintf('%-5s ',names_tmp{p_index(i)})
    for j = 1:i % length(p_index),
        fprintf('%11.4g ',corr_rnd(i,j))
    end
    fprintf('\n')
end
if any(pmat(:,5)==0)
    disp('Note: some parameters have been estimated on log-scale; correlations are also on log-scale for these parameters.')
end
disp('======================================================================================')
disp(' ')
disp(['Number of MCMC samples taken: ',num2str(size(rnd,1))]) % ,', of which accepted in the credible region: ',num2str(size(acc,1))])

% % Checking a for a better optimum is nice but requires calculating a
% % likelihood for every element of rnd (which the slice sampler does, but
% % does not put out). We can calculate it afterwards, but that is time
% % consuming.
% minloglik_max = transfer(parshat,pmat,WRAP); % calculate the best minlog-lik so far ...
% if min(minloglik_rnd)<minloglik_max
%     ind_better = find(minloglik_rnd==min(minloglik_rnd),1,'first');
%     disp(' ')
%     disp(['Warning: a better optimum was found during the MCMC run: ',num2str(min(minloglik_rnd))])
%     for parnr = 1:length(p_index) % run through all fitted parameters
%         if pmat(p_index(parnr),5) == 0 % if parameter was fitted on log-scale, the value in rnd is also on log scale
%             fprintf('parameter %-6s %10.5g \n',names_tmp{p_index(parnr)},10^(rnd(ind_better,parnr)))
%         else
%             fprintf('parameter %-6s %10.5g \n',names_tmp{p_index(parnr)},rnd(ind_better,parnr))
%         end
%     end
%     disp(' ')
% end

disp(['Time required: ' secs2hms(toc)])
diary off    % close results.out

%% Sub-function testslice to test the properties of the sample

function testslice(trace,names)
% Based on code copied from http://www.mbfys.ru.nl/~robvdw/CNP04/LAB_ASSIGMENTS/LAB05_CN05/MATLAB2007b/stats/html/bayesdemo.html

global glo

[nsamples,npars] = size(trace);
scrsz = get(groot,'ScreenSize');
figh = figure;
[~,ft] = make_fig(1,1,1);
lenpl = min([scrsz(4)-100 1500 150*npars]);
widpl = 250*3;
% figh.Position = [300 100 widpl lenpl]; % and give them a default size for mx1
figh.Position = [scrsz(3)-widpl-50 scrsz(4)-lenpl-100 widpl lenpl]; % and give them a default size for 2x2
hold on

ymin = min(trace,[],1);
ymax = max(trace,[],1);
yextr = 0.02*(ymax-ymin);
ybnds = [ymin'-yextr' ymax'+yextr'];

for i = 1:npars
    g = subplot(npars,3,1+3*(i-1));
    plot(trace(:,i)) % plot the random sample in raw form
    xlim([0 length(trace)])
    ylim(ybnds(i,:))
    set(gca,'LineWidth',1,ft.name,ft.ticks)
    p = get(g,'position');
    p([3 4]) = p([3 4]).*[1.0 1.15]; % add some percent to width and height
    if npars == 2 % when plotting just 2 pars, need to move them down a bit
        p(2) = p(2)-0.04;
    end
    set(g, 'position', p);
    hold on
    if i == 1
        title('Raw trace of sample','Interpreter','none',ft.name,ft.title)
    end
    ylabel(names{i},ft.name,ft.label)
    % subplot(length(parshat),1,i),ylabel(names{p_index(i)},'FontSize',10,'FontName','Arial')
    if i<npars
        set(gca,'XTickLabel',[]);
    end
end
xlabel('Number of samples','FontSize',12,'FontName','Arial')

movavg = filter(repmat(1/50,50,1),1,trace); % moving average with window of 50
movavg(1:50,:) = NaN;

for i = 1:npars
    g = subplot(npars,3,2+3*(i-1));
    set(gca,'LineWidth',1,ft.name,ft.ticks)
    p = get(g,'position');
    p([3 4]) = p([3 4]).*[1.05 1.15]; % add some percent to width and height
    if npars == 2 % when plotting just 2 pars, need to move them down a bit
        p(2) = p(2)-0.04;
    end
    set(g, 'position', p);
    hold on
    if i == 1
        title('Moving average of sample','Interpreter','none',ft.name,ft.title)
    end
    plot(movavg(:,i))
    plot([50 50],[min(movavg(:,i)),max(movavg(:,i))],'k:')
    ylim(ybnds(i,:))
    xlim([0 length(trace)])
    % ylim([min(movavg(:,i)),max(movavg(:,i))])
    % ylabel(names{i},'FontSize',10,'FontName','Arial')
    if i<npars
        set(gca,'XTickLabel',[]);
    end
end
xlabel('Number of samples',ft.name,ft.label)

F    =  fft(detrend(trace,'constant'));
F    =  F .* conj(F);
ACF  =  ifft(F);
ACF  =  ACF(1:21,:);                          % Retain lags up to 20.
% ACF  =  real([ACF(1:21,1) ./ ACF(1,1) ...
%              ACF(1:21,2) ./ ACF(1,2)]);       % Normalize.

ACF_tmp = zeros(size(ACF));
for i = 1:npars
    ACF_tmp(:,i) = ACF(1:21,i) ./ ACF(1,i); % normalise
end
ACF = real(ACF_tmp);
    
bounds  =  sqrt(1/nsamples) * [2 ; -2];       % 95% CI for iid normal

for i = 1:npars
    g = subplot(npars,3,3+3*(i-1));
    set(gca,'LineWidth',1,ft.name,ft.ticks)
    p = get(g,'position');
    p([3 4]) = p([3 4]).*[1.15 1.15]; % add some percent to width and height
    if npars == 2 % when plotting just 2 pars, need to move them down a bit
        p(2) = p(2)-0.04;
    end
    set(g, 'position', p);
    hold on
    if i == 1
        title('Autocorrelation in sample','Interpreter','none',ft.name,ft.title)
    end
    lineHandles  =  stem(0:20, ACF(:,i) , 'filled' , 'r-o');
    set(lineHandles , 'MarkerSize' , 4)
    grid('on')
    
    % ylabel(names{i},'FontSize',10,'FontName','Arial')
    hold('on')
    plot([0.5 0.5 ; 20 20] , [bounds([1 1]) bounds([2 2])] , '-b');
    plot([0 20] , [0 0] , '-k');
    hold('off')
    a  =  axis;
    axis([a(1:3) 1]);
    if i<npars
        set(gca,'XTickLabel',[]);
    end
end
xlabel('Lag',ft.name,ft.label)

if glo.saveplt > 0 % if we want to save the plot
    filenm = glo.basenm;
    savenm = ['inspect_slice_',filenm];%
    save_plot(figh,savenm);
end  
