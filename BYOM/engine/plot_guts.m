function plot_guts(par,opt_conf,opt_guts)

% Usage: plot_guts(par,opt_conf,opt_plot)
%
% THIS FUNCTION IS REPLACED BY THE MORE GENERAL plot_tktd FUNCTION.
% HOWEVER, IT IS KEPT IN THE ENGINE FOR NOW AS IT HAS SOME ADDITIONAL
% FUNCTIONALITY NOT (YET) INCLUDED IN THAT NEW FUNCTION. HOWEVER, I HAVE
% NOT TESTED IT WITH THE BYOM ENGINE v.6.0.
% 
% This function provides a series of additional plots, specifically for
% GUTS purposes (analysis of survival data over time). The standard
% plotting routine always plots survival versus time, whereas other
% plotting options may be more useful sometimes: plotting survival versus
% concentration (classic dose-response curves), or plotting observed and
% expected deaths in each observation interval (which is more closely
% representing how the optimisation criterion operates). Furthermore,
% instead of plotting all data in a single plot, this script provides
% multiplots (either different plots for different treatments, or different
% plots for different time points).
%
% Leave opt_conf empty to plot without CIs.
%
% Author     : Tjalling Jager
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html
%
% LIMITATIONS
% - Does not work with multiple data sets yet (only uses first one).
% - Sampling error in dose-response only works properly when number of
% animals is the same for all treatments (otherwise an average is used).
% - Dose-response curves will not work properly if your X0mat entries are
% different for each treatment (first entry is always used).
% - Watch out when using more complex GUTS models, e.g., with time-varying
% exposure, extended TK models, etc. Unexpected results may occur.
% - Several plots will be inaccurate if there are NaNs in the data set.
% - Check results carefully when the exposure concentration is not constant
% over time. The survival-time and deaths-time plots should be OK.

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global DATA W X0mat % make the data set and initial states global variables
global glo % glo2     % allow for global parameters in structure glo

% extrac values from globals
t       = glo.t;
ylab    = glo.ylab;
xlab    = glo.xlab;
leglab1 = glo.leglab1;
leglab2 = glo.leglab2;
% ttot    = glo2.ttot;
% ctot    = glo2.ctot;

% Options for plotting
incl_timeresp  = opt_guts.timeresp; % make a multiplot with survival versus time for each treatment
incl_doseresp  = opt_guts.doseresp; % make a multiplot with survival versus concentration for each time point
incl_timedeath = opt_guts.timedeath; % make a multiplot with deaths-per-interval versus time for each treatment
incl_single_dr = opt_guts.single_dr; % make a single plot with survival versus concentration for each time point
incl_single_td = opt_guts.single_td; % make a single plot with expected and observed deaths in each time interval

bw = opt_guts.bw; % plot the single dose-response plot in black and white
cn = opt_guts.cn; % use connecting lines between data and model curve when opt_guts.bw=1

if isfield(glo,'int_scen') % if it exists, we have defined time-varying concentrations
    if incl_doseresp == 1
        disp(' ')
        disp('No dose-response will be plotted as time-varying exposure treatments are defined.')
    end
    incl_doseresp = 0; % don't do a dose response
end

scriptnm = glo.basenm;

incl_conf   = 0; % do NOT include CIs on plots
incl_samerr = 0; % do NOT include sampling error bounds
if ~isempty(opt_conf) % if left empty, you want to skip all CIs
    incl_conf = 1; % include CIs on plots (but there must be a sample saved)
    if opt_conf.samerr == 1 % if you calculated sampling errors before, I assume you want them here too ...
        incl_samerr = 1; % include sampling error bounds
    end
end 
if incl_conf == 0, type_conf = -1; end

no_c = 50; % number of concentrations in dose-response curves

supp_warn = 0; % don't suppress warning, unless we already caught one

%% Extract survival data, turn it into a probability, and modify for any missing animals
plotdata = DATA{1,glo.locS}; % extract survival data
a        = plotdata(2:end,2:end); % extract data set
w        = W{1,glo.locS}; % extract weight factors

% print warnings, if needed
warning('off','backtrace') % no need to display where the warning is generated
if size(DATA,1)>1
    warning('plot_guts: there are multiple data sets, only first one is plotted!')
end
if any(isnan(a(:)))
    warning('plot_guts: your survival data contains NaNs. Interpret plots carefully!')
end
if sum(w(:))>0
    warning('plot_guts: your survival data contains missing animals. Interpret plots carefully!')
end
warning('on','backtrace')

% code below is from my DEBtoxM code in plot_all. This is a bit
% awkward, but is needed when there are animals missing or
% removed during the experiment. What is plotted is the
% corrected data!
if sum(w(:)) == 0 % if there are NO missing data ...
    a = a ./ (ones(size(a,1),1)*a(1,:)); % turn them into survival probabilities
else % if there are missing data ...
    if ~any(isnan(a(:))) % does the dataset NOT contain nans?
        divideby = (a(1:end-1,:)-w(1:end-1,:));
        divideby(divideby==0) = 1e-6; % to avoid division by zero!
        a = cumprod([ones(1,size(a,2)); a(2:end,:)./divideby]);
    else % it does contain NaNs, which makes life miserable
        w2 = w;
        w2(isnan(a)) = NaN; % also make corresponding points in w2 a NaN
        indnan = cumsum(sum(isfinite(w2),1)); % linear indices to last datapoint for each conc, excl NaNs
        indnan(end) = []; % remove last one (does not overlap with next conc)
        w2 = w2(:); % make them into a column vector
        Z2 = a(:);
        w2(isnan(w2)) = []; % and remove the NaNs
        Z2(isnan(Z2)) = [];
        divideby = (Z2(1:end-1)-w2(1:end-1));
        divideby(divideby==0) = 1e-6; % to avoid division by zero!
        Z3 = Z2(2:end)./divideby; % corrected survival fraction
        Z3(indnan) = []; % remove the overlapping data with next concentration
        Z4 = nan(size(a)); % prepare a new matrix with NaNs
        Z4(1,:) = []; % but take out the first row (will become ones)
        Z4(isfinite(a(2:end,:))) = Z3; % put the corrected fraction back in a matrix at the correct positions
        Z4 = [ones(1,size(a,2));Z4]; % and add a row with ones on top
        Z4(isnan(a))=1; % replace NaNs by one to be able to do cumprod
        Z4 = cumprod(Z4); % cumulative product to get overall survival fraction
        Z4(isnan(a)) = NaN; % put the NaNs back in their position
        a = Z4; % this is the corrected data to plot!
    end
end
n_s = plotdata(2:end,2:end); % remember numbers of survivors as well
plotdata(2:end,2:end) = a; % and put them back
t_d = plotdata(2:end,1); % time vector in the data set
c_d = plotdata(1,2:end); % conc vector in the data set
n_S = a; % survival probability

%% Standard plots of survival versus time
%
% A multiplot is made with survival vs. time for each treatment separately.

if incl_timeresp == 1
    
    c_du = unique(c_d); % unique values of concentrations in data set
    
    if incl_conf ~= 0 || incl_samerr ~= 0 % if we need bounds on curve, invoke calc_conf
        
        opt_conf_rem = opt_conf; % remember options for CIs, so we can put them back later
        X0mat_rem    = X0mat;
        
        % could be that X0mat does not match data set, and we will strictly
        % follow the data set!
        X0mat = zeros(size(X0mat,1),length(c_du));
        X0mat(1,:) = c_du;
        for jc = 1:length(c_du)
            ind_cdu = X0mat(1,:)==c_du(jc);
            X0mat(2:end,jc) = X0mat_rem(2:end,ind_cdu);
        end
        
        opt_conf.sens = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
        
        % use calc)conf to calculate the bounds on the dose-response curve
        % out_conf = calc_conf(par,opt_conf,opt_plot.repls);
        out_conf = calc_conf(par,opt_conf,1,supp_warn); % TEMP ALWAYS PLOT REPLICATES
        supp_warn = 1; % suppress all following warnings from calc_conf
        
        % extract the bounds from out_conf
        Xlo  = out_conf{1};
        Xhi  = out_conf{2};
        XloS = out_conf{3};
        XhiS = out_conf{4};
        type_conf = out_conf{8};
        opt_conf  = opt_conf_rem; % put back the original version
        X0mat     = X0mat_rem;
    end
    
    % initialise the matrices that will collect the output (time x concentrations)
    % we collect only one scenario! (glo.locS)
    X_dr = zeros(length(t),length(c_du)); % initialise state variable output with zeros
    
    for jc = 1:length(c_du)
        ind_cdu = X0mat(1,:)==c_du(jc);
        Xout = call_deri(t,par,X0mat(:,ind_cdu),glo); % use call_deri.m to provide the output for one scenario
        X_dr(:,jc) = Xout(:,glo.locS); % collect survival into structure X_dr
    end
    
    coll_error = []; % initialise to make an overall predicted-observed plot
    
    % Calculate size of multiplot
    n = ceil(sqrt(length(c_du)));
    m = ceil(length(c_du)/n);
    [~,ft] = make_fig(m,n); % create a figure window of correct size
    hold on
    for jc = 1:length(c_du)
        h_pl = subplot(m,n,jc);
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        
        if m>1 % only shrink white space when there are more than 2 rows
            p = get(h_pl,'position');
            p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
            set(h_pl, 'position', p);
        end
        
        if jc>(n*(m-1)) % only put xlabel on bottom row
            xlabel(xlab,ft.name,ft.label)
        else
            set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
        end
        if jc==1 || (jc-1)/n == floor((jc-1)/n) % only put y labels on first column
            ylabel(ylab{glo.locS},ft.name,ft.label)
        else
            set(h_ax,'YTickLabel',[]); % remove tick labels on x-axis
        end
        
        if isfield(glo,'LabelTable') && ismember(c_du(jc),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
            Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == c_du(jc)}; % look up the label belonging to the j-th scenario
            title(Ltmp,ft.name,ft.title)
        else
            title([leglab1,num2str(c_du(jc)),' ',leglab2],ft.name,ft.title)
        end
         
        xlim([0 max(t)])
        ylim([0 1.05])
        plot(t,X_dr(:,jc),'k-','LineWidth',1)
        
        ind_cd = find(c_d == c_du(jc)); % where is this concentration in the data?
        plot(t_d,n_S(:,ind_cd),'ko','MarkerFaceColor','k')
        
        % add values for making an overall predicted-observed plot
        % I assume that the first entry is t=0, so that one is excluded
        for jcr = 1:length(ind_cd)
            coll_error = [coll_error;[n_S(2:end,ind_cd(jcr)) interp1(t(2:end),X_dr(2:end,jc),t_d(2:end))]];
        end
        
        if incl_conf ~= 0 && exist('Xlo','var') % do we have max and min model lines from the slice sampling?
            plot(t,Xlo{glo.locS}(:,jc),'k:','LineWidth',1) % plot each confidence line as broken line
            plot(t,Xhi{glo.locS}(:,jc),'k:','LineWidth',1) % plot each confidence line as broken line
        end
        if incl_samerr ~= 0 && exist('XloS','var') && numel(XloS)>1 % do we have max and min model lines with sampling error from the slice sampling?
            plot(t,XloS(:,jc),'k--','LineWidth',0.5) % plot each confidence line as broken line
            plot(t,XhiS(:,jc),'k--','LineWidth',0.5) % plot each confidence line as broken line
        end
    end
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    
    switch type_conf
        case 0
            h_txt = text(0.5, 1,'Time-response plot','HorizontalAlignment','center','VerticalAlignment', 'top');
        case -1
            h_txt = text(0.5, 1,'Time-response plot','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 1
            h_txt = text(0.5, 1,'Time-response plot, 95% CI from posterior (slicesample)','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 2
            h_txt = text(0.5, 1,'Time-response plot, 95% CI from lik.-region (likregion)','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 3
            h_txt = text(0.5, 1,'Time-response plot, 95% CI from lik.-region (parspace)','HorizontalAlignment','center','VerticalAlignment', 'top');
    end
    set(h_txt,ft.name,ft.text); % use standard formatting for this header
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['timeresp_',scriptnm];%
        save_plot(gcf,savenm,h_txt);
    end
    
    % Make a plot for observed versus predicted of all data points,
    % excluding t=0
    [~,ft] = make_fig(1,1);
    h_ax = gca;
    set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
    hold on
    plot(coll_error(:,1),coll_error(:,2),'k*')
    plot([0 1],[0 1],'k-')
    plot([0 0.9],[0.1 1],'k:')
    plot([0.1 1],[0 0.9],'k:')
    h_leg = legend('data-model','1:1 line','10% deviation','Location','southeast');
    set(h_leg,ft.name,ft.legend);
    xlabel('observed survival probability',ft.name,ft.label)
    ylabel('predicted survival probability',ft.name,ft.label)
    
    error_surv = abs(coll_error(:,2)-coll_error(:,1)); % absolute difference between observed and predicted
%     disp('===================================================================')
%     disp(['mean absolute difference between observed and predicted: ',num2str(mean(error_surv))])
%     disp(['max. absolute difference between observed and predicted: ',num2str(max(error_surv))])
%     disp(['number of differences larger than 0.1: ',num2str(sum(error_surv>0.1))])
%     disp('===================================================================')
    
    % create a text box with summary in the plot
    dim = [0.15 0.6 0.3 0.3];
    str = {['mean diff.: ',num2str(mean(error_surv))],['max. diff.: ',num2str(max(error_surv))],['number >0.1: ',num2str(sum(error_surv>0.1)),' / ',num2str(size(error_surv,1))]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on',ft.name,ft.annot);
    
end

%% Deaths per interval
%
% Create plots with expected and observed numbers of deaths in each
% observation interval. This is closest to how the optimisation routine
% interprets the goodness of fit. First, an overall plot is made, and next
% sub-plots for each treatment separately.

if incl_timedeath == 1 
    % first check if there are replicates, and if they differ in number of
    % individuals between treatments
    c_du = unique(c_d); % unique values of concentrations in data set
    n_su = zeros(1,length(c_du));
    for jc = 1:length(c_du)
        a = n_s(1,c_d==c_du(jc)); % number of individuals in each replicate
        if sum(diff(a)) ~= 0 % then the number of individuals differs between replicates)
            incl_timedeath = 0;
        end
        n_su(jc) = a(1); % number of individuals over replicates in this treatment
    end
    if incl_timedeath == 0
        warning('off','backtrace') % no need to display where the warning is generated
        warning('I will skip plot of deaths for each interval as the initial number of individuals differs between replicates.')
        warning('on','backtrace')
    end
end

if incl_timedeath == 1
    
    n_d = -diff(n_s,1,1); % number of deaths in each interval (excl. after end of test)
    w = W{1,glo.locS}; % also collect weights ...
    if sum(w(:)) > 0
        error('This calculation is not yet accounting for missing/removed animals.')
    end
    
    % initialise the matrices that will collect the output (time x concentrations)
    % we collect only one scenario! (glo.locS)
    X_dr = zeros(length(t_d),length(c_du)); % initialise state variable output with zeros
    
    for jc = 1:length(c_du)
        ind_cd = X0mat(1,:)==c_du(jc); % find correct column in X0mat
        Xout = call_deri(t_d,par,X0mat(:,ind_cd),glo); % use call_deri.m to provide the output for one scenario
        X_dr(:,jc) = Xout(:,glo.locS); % collect survival probabilities into structure X_dr
    end
    
    p_d = -diff(X_dr,1,1); % modelled probability of deaths in each interval (excl. after end of test)
    e_d = p_d .* (ones(length(t_d)-1,1)*n_su); % expected deaths in each interval
    
    if incl_single_td == 1 % want to plot a single expected vs. observed plot?
        
        markers_dr = {'s','^','o','d','v','p','>','<'}; % complete set of 8 markers
        
        h_leg = [];
        L = [];
        [~,ft] = make_fig(1,2);
        subplot(1,2,1)
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        xlabel('observed no. deaths',ft.name,ft.label)
        ylabel('expected no. deaths',ft.name,ft.label)
        hold on
        for jc = 1:length(c_du)
            ind_cd = find(c_d == c_du(jc));
            for j_u = 1:length(ind_cd)
                if jc <= 8 % for BW plots, handle h is used for the legend, h_top for putting symbols on top
                    htmp = plot(n_d(:,ind_cd(j_u)),e_d(:,jc),['k',markers_dr{jc}],'MarkerFaceColor','w');
                elseif jc <= 16
                    htmp = plot(n_d(:,ind_cd(j_u)),e_d(:,jc),['k',markers_dr{jc-8}],'MarkerFaceColor','k');
                elseif jc <= 24
                    htmp = plot(n_d(:,ind_cd(j_u)),e_d(:,jc),['k',markers_dr{jc-16}],'MarkerFaceColor','y');
                end
            end
            h_leg(jc) = htmp(1); % needed for replicates, to get 1 handle number for the legend
            
            if isfield(glo,'LabelTable') && ismember(c_du(jc),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
                L{jc} = glo.LabelTable.Label{glo.LabelTable.Scenario == c_du(jc)}; % look up the label belonging to the j-th scenario
            else
                L{jc} = [leglab1,num2str(c_du(jc)),' ',leglab2];
            end
            
        end
        
        pl_m = max(max(n_d(:)),max(e_d(:)));
        % corrcoef(n_d(:),e_d(:)); % as crude goodness of fit measure?
        plot([0 pl_m],[0 pl_m],'k:') % plot 1:1 line in there
        h_mv = legend(h_leg,L); % remember handle as we will move legend later
        set(h_mv,ft.name,ft.legend);
        
        % Move the legend to the free space
        h_annot = subplot(1,2,2); % make an extra subplot, and remember the handle
        h1 = gca; % remember the current axis
        dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
        h_mv_tmp = h_mv; % take legend
        dimleg = h_mv_tmp.Position; % position of LAST legend made
        % move it to the new subplot and set it to top-left of panel
        h_mv_tmp.Position = [dim(1) dim(2)+dim(4)-dimleg(4) dimleg(3:4)];
        delete(h1) % throw away the axes of the plot, the legend is enough
        
    end
    
    % separate plots per treatment, but first calculate CIs
    if incl_conf ~= 0 || incl_samerr ~= 0 % if we need bounds on curve, call calc_conf
        
        opt_conf_rem = opt_conf; % remember options for CIs, so we can put them back later
        X0mat_rem    = X0mat;
        t_rem        = glo.t;
        
        glo.t     = t_d; % use only the times we need for the dose-response plot
        % X0mat = [c_d;X0mat(2:end,1) * ones(1,length(c_d))]; % create a new X0mat from c_d
        
        % could be that X0mat does not match data set, and we will strictly
        % follow the data set!
        X0mat = zeros(size(X0mat,1),length(c_du));
        X0mat(1,:) = c_du;
        for jc = 1:length(c_du)
            ind_cdu = X0mat(1,:)==c_du(jc);
            X0mat(2:end,jc) = X0mat_rem(2:end,ind_cdu);
        end
               
        opt_conf.sens     = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
        
        % use calc_conf to calculate the bounds on the dose-response curve
        % out_conf = calc_conf(par,opt_conf,opt_plot.repls,supp_warn);
        out_conf = calc_conf(par,opt_conf,1,supp_warn); % TEMP ALWAYS PLOT REPLICATES
        supp_warn = 1; % suppress all following warnings from calc_conf
        
        % extract the bounds from out_conf
        X2 = out_conf{7};
        XloD = out_conf{5};
        XhiD = out_conf{6};
        type_conf = out_conf{8};
        
        opt_conf = opt_conf_rem; % put back the original version
        X0mat    = X0mat_rem;
        glo.t    = t_rem;
        
        X_CI = X2{glo.locS};
        pd_CI = zeros(length(t_d)-1,length(c_du),size(X_CI,3));
        for k = 1:size(X_CI,3) % run through all samples
            pd_CI(:,:,k) = -diff(X_CI(:,:,k),1,1); % probability of deaths in each interval (incl. after end of test)
        end
        
        if opt_conf.type == 1 % slice sampler
            pd_min = prctile(pd_CI,2.5,3);
            pd_max = prctile(pd_CI,97.5,3);
        else % lik region or 95%-max posterior or grid-region
            pd_min = min(pd_CI,[],3);
            pd_max = max(pd_CI,[],3);
        end
        
        ed_min = pd_min .* (ones(length(t_d)-1,1)*n_su); % expected deaths in each interval
        ed_max = pd_max .* (ones(length(t_d)-1,1)*n_su); % expected deaths in each interval
    end
    
    % Calculate size of multiplot
    n = ceil(sqrt(length(c_du)));
    m = ceil(length(c_du)/n);
    
    % Make a new time vector for the symbols. They are plotted halfway into the
    % observation interval, but expected and observed are slightly shifted so
    % they don't overlap completely.
    t_plot_e = mean([t_d(1:end-1) t_d(2:end)],2)-mean(diff(t_d))/20;
    t_plot_d = mean([t_d(1:end-1) t_d(2:end)],2)+mean(diff(t_d))/20;
    
    if incl_samerr ~= 0 && exist('XloD','var') && numel(XloD)>1
        XloD(end,:) = []; % remove last interval (we're not plotting deaths after the test)
        XhiD(end,:) = []; % remove last interval (we're not plotting deaths after the test)
    end
    
    % find the maximum for the y-axis, as the max of all plotting events
    max_x = max(max(n_d(:)),max(e_d(:)));
    if incl_conf ~= 0 && exist('Xlo','var')
        max_x = max(max_x,max(ed_max(:))); 
    end
    if incl_samerr ~= 0 && exist('XloD','var') && numel(XloD)>1
        max_x = max(max_x,max(XhiD(:)));
    end
    
    [~,ft] = make_fig(m,n); % create a figure window of correct size
    hold on
    for jc = 1:length(c_du)
        h_pl = subplot(m,n,jc);
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        
        if m>1 % only shrink white space when there are more than 2 rows
            p = get(h_pl,'position');
            p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
            set(h_pl, 'position', p);
        end
        
        ind_cd = c_d == c_du(jc); % where is this concentration in the data?
        plot(t_plot_d,n_d(:,ind_cd),'ko','MarkerFaceColor','k') % plot observed deaths
        
        if incl_conf ~= 0 && exist('X2','var')% plot model error as error bar on expected deaths
            for jt = 1:length(t_plot_e)
                plot([t_plot_e(jt) t_plot_e(jt)],[ed_min(jt,jc) ed_max(jt,jc)],'k-','LineWidth',1)
                % plot([t_plot_e(jt) t_plot_e(jt)],[XloD(jt,jc) XhiD(jt,jc)],'k:')
            end
        end
        % plot the sampling error as a line over time
        if incl_samerr ~= 0 && exist('XloD','var') && numel(XloD)>1
            plot(t_plot_e,XloD(:,jc),'k--','LineWidth',1)
            plot(t_plot_e,XhiD(:,jc),'k--','LineWidth',1)
        end
        plot(t_plot_e,e_d(:,jc),'ks:','MarkerFaceColor','w','LineWidth',1) % plot expected value
        
        if isfield(glo,'LabelTable') && ismember(c_du(jc),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
            Ltmp = glo.LabelTable.Label{glo.LabelTable.Scenario == c_du(jc)}; % look up the label belonging to the j-th scenario
            title(Ltmp,ft.name,ft.title)
        else
            title([leglab1,num2str(c_du(jc)),' ',leglab2],ft.name,ft.title)
        end
        
        xlim([0 max(t_d)])
        % ylim([0 max(1,max([n_d(:,jc);e_d(:,jc)]))])
        ylim([0 max_x])
        if jc>(n*(m-1)) % only put xlabel on bottom row
            xlabel('time intervals',ft.name,ft.label)
        else
            set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
        end
        if jc==1 || (jc-1)/n == floor((jc-1)/n) % only put y labels on first column
            ylabel('deaths per interval',ft.name,ft.label)
        else
            set(h_ax,'YTickLabel',[]); % remove tick labels on x-axis
        end
    end
    
    % The goodness of fit tests are commented out. The first test now also
    % works for replicated treatments. However, they seem to be of limited
    % use in judging a fit. 
    
%     % ==============================================================
%     % Test for goodness-of-fit from Monte-Carlo simulations
%     % This does not account for missing/removed animals yet, but in those
%     % cases, an error will be produced above anyway. This is based on
%     % Albert et al, 2012 (DOI 10.1016/j.jspi.2011.07.014).
%     fit_coll = zeros(10000,1); % prepare to collect chi-square criteria for the simulations
%     fit_stat = 0; % this will be the goodness-of-fit statistic for the best fit of the real data
%     for jc = 1:length(c_du) % run through unique treatments
%         ind_cd = find(c_d == c_du(jc)); % find if there are replicates for this one
%         a = sum(n_s(1,ind_cd)); % sum of initial individuals in each replicate
%         b = sum(n_d(:,ind_cd),2); % number of deaths over time, summed over replicates for this treatment
%         c = a * (e_d(:,jc)/n_su(jc)); % convert expected deaths (per replicate) back to probability, back to sum of replicates
%         
%         stat_obs  = [b;a - sum(b)]; % observed deaths (add last time interval, after test)
%         stat_pred = [c ; a - sum(c)]; % predicted deaths (add last time interval, after test)
%         stat_pred = max(1e-20,stat_pred); % avoid zero, which would lead to NaNs for fitstat
%         stat_pd   = [p_d(:,jc) ; 1 - sum(p_d(:,jc))]; % probability of death
%         
%         % add the contributions for this treatment to the fit statistic for the actual data set
%         fit_stat = fit_stat + sum(((stat_obs-stat_pred).^2) ./ stat_pred);
%         
%         % run the simulations for this particular treatment
%         for j_stat = 1:length(fit_coll)
%             % simulate an outcome in terms of numbers of deaths, from the death probabilities
%             deaths_rnd = mnrnd(a,stat_pd);
%             % add the contributions for this treatment to the fit statistic
%             fit_coll(j_stat) = fit_coll(j_stat) + sum(((deaths_rnd'-stat_pred).^2) ./ stat_pred);
%         end
%     end
%     
%     % what percentage of samples is worse than the observed deaths? If that's less than 0.05 that's bad.
%     stat_final = sum(fit_coll>fit_stat)/length(fit_coll);
%     % how do degrees of freedom come in? they don't.
%     if isnan(fit_stat) || any(isnan(fit_coll))
%         warning('off','backtrace') % no need to display where the warning is generated
%         warning('plot_guts: cannot test goodness-of-fit due to NaNs.')
%         warning('on','backtrace'), disp(' ')
%     else
%         disp('============================================================================')
%         % disp(['Goodness-of-fit, chi-square criterium: ', num2str(fit_stat)])
%         disp('Test for goodness-of-fit (Pearson chi-square with MC simulated distribution).')
%         disp(['  Observed fit is better than this fraction of samples (p-value): ', num2str(stat_final)])
%         if stat_final < 0.05
%             disp('  Fit rejected at a=0.05; the data are not consistent with the calibrated model.')
%         else
%             disp('  Fit is not rejected at a=0.05.')
%         end
%         disp('============================================================================')
%     end
    
    % ==============================================================
    
%     % ==============================================================
%     % another goodness-of-fit criterion: the deviance as used by Bedaux and
%     % Kooijman.
%     % THIS DOES NOT WORK WHEN THERE ARE REPLICATES OR NaNs!
%     n_tmp = ones(length(t_d),1)*n_s(1,:); % matrix with initial nr animals per treatement, copied to all time points
%     n_dt = [-diff(n_s,1,1);n_s(1,:)-sum(n_d)]; % number of deaths in each interval (INCL. after end of test)
%     prob_sup = n_dt./n_tmp; % deasth probabilities from the data set
%     prob_sup = max(prob_sup,1e-50); % Otherwise we get problems taking the logarithm
%     loglik_sup = (n_dt(:))'*log(prob_sup(:));  % and calculate the log likelihood
%     prob_mod = [p_d ; 1-sum(p_d)]; % add a row with predicted death probability in interval after the test
%     loglik_mod = (n_dt(:))'*log(prob_mod(:));  % and calculate the log likelihood
%     pmat = packunpack(1,par,0);  % transform structure into a regular matrix
%     loglikrat = -2 * (loglik_mod-loglik_sup); % likelihood ratio criterion
%     GoodFit = 1-chi2cdf(loglikrat,numel(prob_sup)-length(c_d)-sum(pmat(:,2)==1));
%     % df is difference in estimated parameters; for the saturated model,
%     % that is the number of elements in prob_sup with the last row (as the
%     % probabilities need to sum to 1 for each column).
%     
%     disp('============================================================================')
%     disp('Test for goodness-of-fit using the deviance.')
%     disp(['  approximate p-value from chi-square distribution: ', num2str(GoodFit)])
%     if GoodFit < 0.05
%         disp('  Fit rejected at a=0.05.')
%     else
%         disp('  Fit is not rejected at a=0.05.')
%     end
%     disp('============================================================================')
%     
% %     fit_coll = zeros(10000,1); % prepare to collect chi-square criteria for the simulations
% %     for jc = 1:length(c_d) % run through treatments
% %         % run the simulations for this particular treatment
% %         for j_stat = 1:length(fit_coll)
% %             % simulate an outcome in terms of numbers of deaths, from the death probabilities
% %             deaths_rnd = mnrnd(n_s(1,jc),prob_mod(:,jc));
% %             % add the contributions for this treatment to the fit statistic
% %             prob_sup = deaths_rnd'./n_tmp(:,jc); % deasth probabilities from simulated data set
% %             prob_sup = max(prob_sup,1e-50); % Otherwise we get problems taking the logarithm
% %             loglik_sup = deaths_rnd*log(prob_sup);  % and calculate the log likelihood
% %             
% %             fit_coll(j_stat) = fit_coll(j_stat) + loglik_sup;
% %         end
% %     end
%     
%     
%     % ==============================================================
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    switch type_conf
        case -1
            h_txt = text(0.5, 1,'Deaths per interval','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 1
            h_txt = text(0.5, 1,'Deaths per interval, 95% CI from posterior (slicesample)','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 2
            h_txt = text(0.5, 1,'Deaths per interval, 95% CI from lik.-region (likregion)','HorizontalAlignment','center','VerticalAlignment', 'top');          
        case 3
            h_txt = text(0.5, 1,'Deaths per interval, 95% CI from lik.-region (parspace)','HorizontalAlignment','center','VerticalAlignment', 'top');          
    end
    set(h_txt,ft.name,ft.text); % use standard formatting for this header
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['deaths_',scriptnm];%
        save_plot(gcf,savenm,h_txt);
    end
    
end

%% Dose-response curves
%
% Create plots with survival vs. concentration. First a single plot is made
% with all requested time points. Next, a multiplot with subplots for each
% time point in the data set.

if incl_doseresp == 1
    
    cd_rem = c_d; % remember the original concentration range, as we will modify it
    
    % Use a log-range with as minimum the lowest non-zero concentration divided
    % by the average step in the data set. This is needed to get nice plots on
    % log scale. The control (if it is zero) is plotted separately, which its
    % own tiny range (contr_range).
    c_du = unique(c_d); % only unique values, and sort
    c_step = mean(diff(log10(c_du(c_du>0)))); % on log scale
    c_range = logspace(log10(min(c_du(c_du>0)))-c_step,log10(max(c_du)),no_c);
    if min(cd_rem) == 0
        contr_range = [min(c_range)/(10^(1.5*c_step)) min(c_range)/(10^(0.5*c_step)) ];
    else
        contr_range = [NaN NaN];
    end
    c_range = [0 c_range]; % make sure there is always a zero in the range
    
    % initialise the matrices that will collect the output (time x concentrations)
    % we collect only one scenario! (glo.locS)
    X_dr = zeros(length(t_d),length(c_range)); % initialise state variable output with zeros
    
    for jc = 1:length(c_range) % run through large concentration range
        % take first column in X0mat for all concentrations!
        Xout = call_deri(t_d,par,[c_range(jc); X0mat(2:end,1)],glo); % use call_deri.m to provide the output for one scenario
        X_dr(:,jc) = Xout(:,glo.locS); % collect survival into structure X_dr
    end
    
    if incl_conf ~= 0 || incl_samerr ~= 0 % if we need bounds on curve, run calc_conf
        
        if incl_samerr == 1 && sum(n_s(1,:)==mean(n_s(1,:)))~=length(n_s(1,:))
            % then the number of animals differs per treatment
            disp(' '), warning('off','backtrace')
            warning('The number of animals per treatment differs. The bounds including sampling error will be innaccurate (they are now based on the mean number of animals over all replicated).')
            warning('off','backtrace'), disp(' ')
        end
        
        opt_conf_rem = opt_conf; % remember options for CIs, so we can put them back later
        X0mat_rem    = X0mat;
        t_rem        = glo.t;
        
        t     = t_d; % use only the times we need for the dose-response plot
        rem_1 = 0;
        if min(t) ~= min(t_rem)
            t = [min(t_rem) t]; % add time zero in there, this is needed for the sample error!
            rem_1 = 1; % identifier that we included an extra time point
        end
        glo.t = t; % change the global to the limited time vector!
        X0mat = [c_range;X0mat(2:end,1) * ones(1,length(c_range))]; % create a new X0mat from c_range
        
        opt_conf.sens     = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
        
        % use calc)conf to calculate the bounds on the dose-response curve
        % out_conf = calc_conf(par,opt_conf,opt_plot.repls,supp_warn);
        out_conf = calc_conf(par,opt_conf,1,supp_warn,1); % TEMP FOR NOW ALWAYS DO REPLICATES
        % last input is to always use means over all replicates/treatments
        % for numbers of individuals when calculating sampling errors
        supp_warn = 1; % suppress all following warnings from calc_conf
        
        % extract the bounds from out_conf
        Xlo = out_conf{1};
        Xhi = out_conf{2};
        XloS = out_conf{3};
        XhiS = out_conf{4};
        type_conf = out_conf{8};
        
        % Xlo and Xhi are now calculated
        if rem_1 == 1 % need to remove the first time point
            Xlo{glo.locS}(1,:) = [];
            Xhi{glo.locS}(1,:) = [];
            if opt_conf.samerr ~= 0
                XloS(1,:) = [];
                XhiS(1,:) = [];
            end
        end
        
        opt_cont = opt_conf_rem; % put back the original version
        X0mat    = X0mat_rem;
        glo.t    = t_rem;
        t = t_rem;
    end
    
    if incl_single_dr == 1 % want a single plot?
        % first plot a single figure with the times requested in t_d
        plotcol = 'krbgcm'; % 6 colors should be enough ...
        markers_dr = {'s','^','o','d','v','p','>','<'}; % complete set of 8 markers
        
        [~,ft] = make_fig(1,2);
        subplot(1,2,1)
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        xlabel([leglab1,' (',leglab2,')'],ft.name,ft.label)
        ylabel(ylab{glo.locS},ft.name,ft.label)
        
        L = [];
        h_leg = [];
        for jt = 1:length(t_d)
            if bw == 1 % plot in black and white
                plot(c_range,X_dr(jt,:),'k-','LineWidth',1) % plot each model curve with a black line
                plot(contr_range,[X_dr(jt,1) X_dr(jt,1)],'k-','LineWidth',1) % plot control range
            else
                if jt <= 6
                    htmp = plot(c_range,X_dr(jt,:),[plotcol(jt),'-'],'LineWidth',2); % plot each model curve with a different color
                    plot(contr_range,[X_dr(jt,1) X_dr(jt,1)],[plotcol(jt),'-'],'LineWidth',2) % plot control range
                elseif jt <= 12 % if we have more than 6 treatments, re-use colours and make broken line
                    htmp = plot(c_range,X_dr(jt,:),[plotcol(jt-6),'--'],'LineWidth',2); % plot each model curve with a different color
                    plot(contr_range,[X_dr(jt,1) X_dr(jt,1)],[plotcol(jt-6),'-'],'LineWidth',2) % plot control range
                elseif jt <= 18 % if we have more than 12 treatments, re-use colours and make dot-stripe line
                    htmp = plot(c_range,X_dr(jt,:),[plotcol(jt-12),'.-'],'LineWidth',2); % plot each model curve with a different color
                    plot(contr_range,[X_dr(jt,1) X_dr(jt,1)],[plotcol(jt-12),'-'],'LineWidth',2) % plot control range
                else
                    disp('cannot plot more than 18 lines yet ... so modify your script')
                end
                h_leg(jt) = htmp(1);
                L{jt} = [xlab,' ',num2str(t_d(jt))];
            end
        end
        
        if bw == 0
            h_mv = legend(h_leg,L); % remember handle as we will move legend later
        end
        
        % We will use log-scale, so plot the control at geometric mean of control range
        if min(cd_rem) == 0 % if there is a zero in there ...
            c_d(c_d==0) = geomean(contr_range); % use mean of control range (>0)
        end
        
        ylim([0 1.05]); % restrict axis from 0-1
        h_leg = [];
        count_lines = 0; % counter for number of lines plotted (for legend)
        for jt = 1:length(t_d)
            i_td = find(t_d == t_d(jt)); % find this time point in the data set
            if ~isempty(i_td)
                count_lines = count_lines + 1;
                if bw == 1 % plot in black and white
                    if cn == 1 % plot connection between line and point
                        plotcon = interp1(c_range,X_dr(jt,:),c_d); % model prediction at measured data points
                        for con = 1:length(c_d)
                            if ~isnan(n_S(i_td,con)) % only do that when there is a data point
                                plot([c_d(con) c_d(con)],[plotcon(con) n_S(i_td,con)],'k:','LineWidth',1);
                            end
                        end
                    end
                    
                    if jt <= 8 % for BW plots, handle h is used for the legend, h_top for putting symbols on top
                        htmp = plot(c_d,n_S(i_td,:),['k',markers_dr{jt}],'MarkerFaceColor','w');
                    elseif jt <= 16
                        htmp = plot(c_d,n_S(i_td,:),['k',markers_dr{jt-8}],'MarkerFaceColor','k');
                    elseif jt <= 24
                        htmp = plot(c_d,n_S(i_td,:),['k',markers_dr{jt-16}],'MarkerFaceColor','y');
                    end
                    h_leg(count_lines) = htmp(1); % needed for replicates, to get 1 handle number for the legend
                    L{count_lines} = [xlab,' ',num2str(t_d(jt))];
                    
                else % for color plots, handle h_top is used only for putting symbols on top
                    if jt <= 6
                        htmp = plot(c_d,n_S(i_td,:),'ko','MarkerFaceColor',plotcol(jt));
                    elseif jt <= 12 % if we have more than 6 treatments, re-use colours and make squares
                        htmp = plot(c_d,n_S(i_td,:),'ks','MarkerFaceColor',plotcol(jt-6));
                    elseif jt <= 18 % if we have more than 12 treatments, re-use colours and make triangles
                        htmp = plot(c_d,n_S(i_td,:),'k^','MarkerFaceColor',plotcol(jt-12));
                    end
                    h_leg(count_lines) = htmp(1); % needed for replicates, to get 1 handle number for the legend
                end
                
            end
        end
        
        if bw == 1
            h_mv = legend(h_leg,L); % remember handle as we will move legend later
            set(h_mv,ft.name,ft.legend);
        end
        drawnow
        
        if min(cd_rem) == 0
            xlim([min(contr_range) max(c_range)])
            plot([contr_range(2) contr_range(2)],[0 1.02],'b-','LineWidth',0.5)
            plot([c_range(2) c_range(2)],[0 1.02],'b-','LineWidth',0.5)
        else
            xlim([min(c_range(2:end)) max(c_range(2:end))]) % because we added a zero!
        end
        
        h_ax.XScale = 'log'; % put x-axis on log-scale
        uistack(h_leg,'top') % put symbols on top, at the end of all plotting
        
        % Move the legend to the free space in the next subplot
        h_annot = subplot(1,2,2); % make an extra subplot, and remember the handle
        h1 = gca; % remember the current axis
        dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
        h_mv_tmp = h_mv; % take legend
        dimleg = h_mv_tmp.Position; % position of LAST legend made
        % move it to the new subplot and set it to top-left of panel
        h_mv_tmp.Position = [dim(1) dim(2)+dim(4)-dimleg(4) dimleg(3:4)];
        delete(h1) % throw away the axes of the plot, the legend is enough
        drawnow
        
        if glo.saveplt > 0 % if we want to save the plot
            savenm = ['doseresp_single_',scriptnm];%
            save_plot(gcf,savenm);
        end
        
    end
    
    %% Dose-response curves in multiplot
    %
    % Create plots with survival versus concentration; one plot for each time
    % point in the data set.
    
    % Calculate size of the multiplot
    n = ceil(sqrt(length(t_d)-1));
    m = ceil((length(t_d)-1)/n);
    
    [~,ft] = make_fig(m,n); % create a figure window of correct size
    hold on
    for jt = 2:length(t_d)
        h_pl = subplot(m,n,jt-1);
        hold on
        h_ax = gca;
        set(h_ax,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        
        if m>1 % only shrink white space when there are more than 2 rows
            p = get(h_pl,'position');
            p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
            set(h_pl, 'position', p);
        end
        
        if jt-1>(n*(m-1)) % only put xlabel on bottom row
            xlabel([leglab1,' (',leglab2,')'],ft.name,ft.label)
        else
            set(h_ax,'XTickLabel',[]); % remove tick labels on x-axis
        end
        if jt==2 || (jt-2)/n == floor((jt-2)/n) % only put y labels on first column
            ylabel(ylab{glo.locS},ft.name,ft.label)
        else
            set(h_ax,'YTickLabel',[]); % remove tick labels on x-axis
        end
        
        plot(c_range,X_dr(jt,:),'k-','LineWidth',1) % plot each model curve with a black line
        plot(contr_range,[X_dr(jt,1) X_dr(jt,1)],'k-','LineWidth',1) % plot control range
        
        if incl_conf ~= 0 && exist('Xlo','var') % do we have max and min model lines from the slice sampling?
            plot(c_range,Xlo{glo.locS}(jt,:),'k:','LineWidth',1) % plot each confidence line as broken line
            plot(c_range,Xhi{glo.locS}(jt,:),'k:','LineWidth',1) % plot each confidence line as broken line
            if min(cd_rem) == 0
                plot(contr_range,[Xlo{glo.locS}(jt,1) Xlo{glo.locS}(jt,1)],'k:','LineWidth',1) % plot control
                plot(contr_range,[Xhi{glo.locS}(jt,1) Xhi{glo.locS}(jt,1)],'k:','LineWidth',1) % plot control
            end
        end
        if incl_samerr ~= 0 && exist('XloS','var') && numel(XloS)>1 % do we have max and min model lines with sampling error from the slice sampling?
            plot(c_range,XloS(jt,:),'k--','LineWidth',0.5) % plot each confidence line as broken line
            plot(c_range,XhiS(jt,:),'k--','LineWidth',0.5) % plot each confidence line as broken line
            if min(cd_rem) == 0
                plot(contr_range,[XloS(jt,1) XloS(jt,1)],'k--','LineWidth',0.5) % plot control
                plot(contr_range,[XhiS(jt,1) XhiS(jt,1)],'k--','LineWidth',0.5) % plot control
            end
        end
        
        plot(c_d,n_S(jt,:),'ko','MarkerFaceColor','k'); % plot the data
        
        if min(cd_rem) == 0
            plot([contr_range(2) contr_range(2)],[0 1.02],'b-','LineWidth',0.5)
            plot([c_range(2) c_range(2)],[0 1.02],'b-','LineWidth',0.5)
            plot(geomean(contr_range),n_S(jt,cd_rem==0),'ko','MarkerFaceColor','k'); % plot the data
        end
        
        title([xlab,' ',num2str(t_d(jt))],ft.name,ft.title)
        h_ax.XScale = 'log'; % put x-axis on log-scale
        ylim([0 1.05])
        if min(cd_rem) == 0
            xlim([min(contr_range) max(c_range)])
        else
            xlim([min(c_range(2:end)) max(c_range(2:end))]) % because we added a zero!
        end
        
    end
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    switch type_conf
        case -1
            h_txt = text(0.5, 1,'Conc.-response plot','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 1
            h_txt = text(0.5, 1,'Conc.-response plot, 95% CI from posterior (slicesample)','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 2
            h_txt = text(0.5, 1,'Conc.-response plot, 95% CI from lik.-region (likregion)','HorizontalAlignment','center','VerticalAlignment', 'top');
        case 3
            h_txt = text(0.5, 1,'Conc.-response plot, 95% CI from lik.-region (parspace)','HorizontalAlignment','center','VerticalAlignment', 'top');
    end
    set(h_txt,ft.name,ft.text); % use standard formatting for this header
    
    c_d = cd_rem; % return the previous version of the concentration vector
    drawnow
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['doseresp_',scriptnm];%
        save_plot(gcf,savenm,h_txt);
    end
    
end
