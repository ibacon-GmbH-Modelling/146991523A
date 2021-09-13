function [Xall,Xall2x,Xall2y] = calc_and_plot(par_plot,opt_plot,varargin)

% Usage: [Xall,Xall2x,Xall2y] = calc_and_plot(par_plot,opt_plot,varargin)
%
% Fit or not fit, calculate model curves, and plot everything. Optional
% outputs include the model curves in <Xall>, and for the extra data sets
% in <Xall2x> and <Xall2y>. These can be used to make custom plots in your
% main script. <varargin> is used to allow the input of information on
% confidence intervals (output from <calc_conf.m>).
% 
% Author     : Tjalling Jager 
% Date       : March 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global DATA W X0mat % make the data set and initial states global variables
global DATAx Wx     % optional global for extra data set with different x-axis
global glo glo2     % allow for global parameters in structure glo

global h_txt 
global fh_tmp % to save different plots in the main script ...

% read options from structure
zvd_plot = opt_plot.zvd;     % plotting of zero-variate data (if defined)
sub_plot = opt_plot.sub;     % switch for putting state variables as sub-plots into a figure
bw       = opt_plot.bw;      % if set to 1, plots in black and white with different plot symbols
cn       = opt_plot.cn;      % if set to 1, connect model line to points (only for bw=1)
limax    = opt_plot.limax;   % if set to 1, limit axes to the data set for each stage
sho      = opt_plot.sho;     % set to 1 to show all scenarios, 0 to only plot model for scenarios with data
repls    = opt_plot.repls;   % set to 1 to plot replicates, 0 to plot mean responses
outl     = opt_plot.outl;    % set to 1 to identify outliers in the plot (points with weight zero; not for survival)
legsup   = opt_plot.legsup;  % set to 1 to suppress legends on fits 
annot    = opt_plot.annot;   % annotations in sub-plot: text box with parameter estimates or overall legend
simulik  = opt_plot.simulik; % if set to 1, calculate the log-likelihood when not fitting
statsup  = opt_plot.statsup; % vector with states to suppress in plotting fits
notitle  = opt_plot.notitle; % if set to 1, suppress titles
transf   = opt_plot.transf;  % set to 1 to calculate means and SEs including transformations
y_zero   = opt_plot.y_zero;  % set to 1 to force y-axis to start at zero

% extract parameters from the general globals, glo and glo2
savenm1 = glo.basenm; % basis for title and name for saving fitting graphs (when glo.saveplt>0)
t       = glo.t;
n_D     = glo2.n_D;
n_X     = glo2.n_X;
n_X2    = glo2.n_X2;
ttot    = glo2.ttot;
ctot    = glo2.ctot;
names   = glo2.names;
namesp  = glo2.namesp;
namesz  = glo2.namesz;
 
if ~isempty(varargin) && ~isempty(varargin{1}) % then this function is called with info on CI's
    out_conf  = varargin{1};
    Xlo       = out_conf{1};
    Xhi       = out_conf{2};
    XloS      = out_conf{3};
    XhiS      = out_conf{4};
    type_conf = out_conf{8};
end

n_s = size(X0mat,2); % number of scenarios

% prelim_checks is now called from the main script, as it also defines all
% the option sets which need to be in the workspace.

show_X          = 1:n_X; % states to plot
show_X(statsup) = []; % but remove states that are to be suppressed

h_txt = cell(1,n_D); % make sure handle for text titles exists and is empty (to allow suppressing them in save_plot)

if length(show_X) == 1 && annot == 0 % if there is only one state, don't act as if there are sub-plots
    sub_plot = 0; % unless you choose the annotate option
end

%% Derive minloglik if needed

FVAL = NaN; % by default, the minloglik gets a NaN
AIC  = NaN; % and Akaike as well
if ~isfield(par_plot,'tag_fitted') % simulate output, so present parameter values on screen
    pmat_plot = packunpack(1,par_plot,0); % extract the parameters that are used for plotting
    pmat_plot(:,2) = 0; % set all parameters to 'not fitted'
    disp(' ') % also make a short summary on screen
    fprintf('Plots are simulations using the user-provided parameter set. \n');
    fprintf('Parameter values used \n');
    fprintf('=====================================\n');
    nfields = length(names);
    for i = 1:nfields % display all parameters on screen
        fprintf('%-6s %10.4g \n',names{i},pmat_plot(i,1))
    end
    fprintf('=====================================\n');
    if simulik == 1 && numel(ttot) > 1 % ttot will be zero when there are no data! 
        % simulate output, and need to calculate minloglik
        pmat_plot = packunpack(1,par_plot,0); % extract the parameters that are used for plotting
        pmat_plot(:,2) = 0; % set all parameters to 'not fitted'
        pmat_plot(pmat_plot(:,5)==0,1) = log10(pmat_plot(pmat_plot(:,5)==0,1));
        FVAL = transfer([],pmat_plot); % use transfer to obtain min log-likelihood
        fprintf('Minus log-likelihood %1.6g. \n',FVAL)
        if ~isempty(namesp)
            fprintf('  Including contributions from prior distributions. \n')
        end
        if ~isempty(namesz)
            fprintf('  Including contributions from zero-variate data. \n')
        end
        fprintf('=====================================\n');
    end
    fprintf('Filename: %s',glo.basenm); % name of the mydata file
    fprintf('  \n');
    if ~isempty(varargin) % this function is called with CI's!
        disp(' '), warning('off','backtrace')
        warning('Model lines represent simulations from the initial parameter structure par, whereas the confidence intervals result from the settings in the saved parameter set.')
        warning('on','backtrace'), disp(' ')
    end
    AIC = NaN; % Akaiki information criterion only makes sense when things are fitted
else
    % disp(' ') % also make a short summary on screen
    fprintf('Plots result from the optimised parameter values. \n');
    % disp(' ')
    if annot == 1 && sub_plot == 1 && numel(ttot) > 1 % ttot will be zero when there are no data!
        % When plotting parameter values in a subplot, always need to obtain the minloglik
        pmat_plot = packunpack(1,par_plot,0); % extract the parameters that are used for plotting
        pmat_plot(:,2) = 0; % set all parameters to 'not fitted'
        pmat_plot(pmat_plot(:,5)==0,1) = log10(pmat_plot(pmat_plot(:,5)==0,1));
        FVAL = transfer([],pmat_plot); % use transfer to obtain min log-likelihood
        pmat_plot = packunpack(1,par_plot,0); % extract again, to get correct 'fitted' tags
        AIC = 2*sum(pmat_plot(:,2))+2*FVAL; % Akaiki information criterion
    end
    
end
    
% on small screens, Matlab does not allow certain figure sizes, so they are
% scaled down in make_fig; just warn the user here that that will happen
scrsz = get(groot,'ScreenSize');
if scrsz(4) < 870 && strcmp(get(0,'DefaultFigureWindowStyle'),'normal')
    % 870 this is probably the limit; above that Matlab starts to scale automatically
    fprintf('\nDue to your small screen size, all multiplots 2x2 and larger will be scaled down.\n\n')
end

%% Calculate model curves for plotting

% initialise the matrices that will collect the output (time x scenario)
for i = 1:n_X % run through state variables
    Xall{i} = zeros(length(t),n_s); % initialise state variable output with zeros
end
% tt = zeros(1,n_s);         % event times
for i = 1:n_X2 % run through additional data sets
    Xall2x{i} = zeros(length(t),n_s); % initialise additional x-data output with zeros
    Xall2y{i} = zeros(length(t),n_s); % initialise additional y-data output with zeros
end

L = cell(n_s,1); % pre-define L as empty cell array

for j = 1:n_s % run through our scenarios
    
    [Xout,TE,Xout2,zvd] = call_deri(t,par_plot,X0mat(:,j),glo); % use call_deri.m to provide the output for one scenario
    % Note: zvd is used in an eval below.
    
    for i = 1:n_X  % run through state variables
        Xall{i}(:,j) = Xout(:,i); % collect state variable i into structure X
    end
    tt{j} = TE; % remember where the event took place
    
    for i = 1:n_X2 % now loop over the additional data sets
        Xall2x{i}(:,j) = Xout2{i}(:,1); % collect the new x-values
        Xall2y{i}(:,j) = Xout2{i}(:,2); % collect the new y-values
    end
    
    % this creates a cell array for the figure legends from the first row of X0mat
    if isfield(glo,'LabelTable') && ismember(X0mat(1,j),glo.LabelTable.Scenario) % there is a definition table for legends, and this scenario is in it
        L{j} = glo.LabelTable.Label{glo.LabelTable.Scenario == X0mat(1,j)}; % look up the label belonging to the j-th scenario
    else    
        L{j} = [glo.leglab1,num2str(X0mat(1,j)),' ',glo.leglab2];
    end
    
end

% NOTE: when there are multiple data sets per state, it is not possible to
% plot additional model curves for scenarios without data ... There does
% not seem to be a simple way around it (that I can think of), so the
% easiest workaround would be to define some dummy data.
if n_D > 1
    % this is to ensure that when there are more data series, that we know
    % which scenarios are in which data series
    for k = 1:n_D % run through number of data sets per state
        c_all = [];
        for i = 1:n_X % run through state variables
            c_all = [c_all DATA{k,i}(1,2:end)];
        end
        ck{k} = unique(c_all);
    end
end

Le = cell(n_D,n_X); % pre-define as empty cell array

for k = 1:n_D % run through number of data sets per state
    for i = 1:n_X % run through state variables
        if sho == 1 % show all scenario's in X0mat
            if n_D > 1
                % when there are more data series, plot only model curves
                % relevant for THAT series
                [dummy,loci]=ismember(ck{k},X0mat(1,:)); % loci is location of state i in the scenario range
                if sum(dummy) ~= 0
                    loci(loci==0) = []; % remove zeros (where data for certain scenarios are ignored)
                    % this seems to be needed for cases with more data sets per state
                    loc{k,i} = loci;
                else
                    loc{k,i} = [];
                end
            else
                loc{k,i} = 1:size(Xall{i},2);
            end
            
        else % otherwise, only show model results for treatments were there are data
            [~,loci]=ismember(DATA{k,i}(1,2:end),X0mat(1,:)); % loci is location of state i in the scenario range
            loci(loci==0) = []; % remove zeros (where data for certain scenarios are ignored)
            loci = unique(loci); % only keep the unique entries (in case of replicated data)
            % X{i} = X{i}(:,loci); % use loci to extract only the model values for scenarios used in data set i
            loc{k,i} = loci; % put the result in a structure for later use
        end
        if ~isempty(loc{k,i})
            Le{k,i} = L(loc{k,i}); % use loc{i} to extract only legend entries for scenarios used in data set i
        else
            Le{k,i} = [];
        end
    end
end
if n_X2 > 0 % if there are extra data sets, we want to know about that for the colours/markers already
    for k =1:n_X2
        [~,loc_tmp]=ismember(DATAx{k}(1,2:end),X0mat(1,:)); % loci is location of state i in the scenario range
        loc_tmp(loc_tmp==0) = []; % remove zeros (where data for certain scenarios are ignored)
        loc2{k} = loc_tmp;
    end
end

%% Plot model curves with data points

plotcol = 'krbgcm'; % 6 colors should be enough ...
markers_dr = {'s','^','o','d','v','p','>','<'}; % complete set of 8 markers

% calculate size of plot for subplots for each state variable
n_X_plot = length(show_X); % new number of states to plot

n = ceil(sqrt(n_X_plot));
m = ceil(n_X_plot/n);
if annot ~= 0 && sub_plot == 1 % make an extra subplot with the parameter estimates or other annotations
    n = ceil(sqrt(n_X_plot+1));
    m = ceil((n_X_plot+1)/n);
end

locall = [];
for k = 1:n_D % run through number of data sets per state
    for i = 1:n_X % run through state variables
        locall = [locall loc{k,i}];
    end
end
if n_X2 > 0 % if there are extra data sets, we want to know about that for the colours/markers already
    for k =1:n_X2 % run through extra data sets
        locall = [locall loc2{k}];
    end
end
locall = unique(locall);

h_leg = []; % clear handles for legends
flag_meanskip = 0; % flag if there is a mean skipped because there are NaNs in some replicates

for k = 1:n_D % run through number of data sets per state
    
     if sub_plot == 1
         [fh_tmp{k},ft] = make_fig(m,n); % create a figure window of correct size
         % I remember them in a global so I can save the various plots in
         % the main script!
     end
    
    for i_tmp = 1:n_X_plot % run through state variables
        i = show_X(i_tmp); % this allows to suppress certain plots!
        
        plotpts = [];
        h     = []; % clear handle for plotting for plot symbols (BW)
        h_top = []; % handles for putting symbols on top
        h_out = []; % handles for putting outlier symbols on top
        
        max_data = 0; % clear maximum for data (to limit axis bounds if needed)
        min_data = +inf; % clear minimum for data (to limit axis bounds if needed)
        
        if sub_plot == 1 % make a subplot for state i
            subplot(m,n,i_tmp)
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        else
            [fig_close,ft] = make_fig(1,1);  % create a figure window of correct size
            fh_tmp{k,i_tmp} = fig_close; % put it in the global for potential later use
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            if notitle == 0
                if ~exist('Xlo','var')
                    if n_D == 1 % add a title to the graph with the filename in it
                        title(['Plotted from: ',savenm1, ' (',date,')'],'Interpreter','none',ft.name,ft.title)
                    else
                        title(['Plotted from: ',savenm1, ' (',date,'), series ',num2str(k),' of ',num2str(n_D)],'Interpreter','none',ft.name,ft.title)
                    end
                else
                    switch type_conf
                        case 1
                            title('Intervals: 95% of predictions from posterior',ft.name,ft.title)
                        case 2
                            title('Intervals: prediction from likelihood-region (likreg)',ft.name,ft.title)
                        case 3
                            title('Intervals: prediction from likelihood-region (parspac)',ft.name,ft.title)
                    end
                end
            end
        end
        hold on
       
        % add labels to the axes
        xlabel(glo.xlab,ft.name,ft.label)
        ylabel(glo.ylab{i},ft.name,ft.label)
        
        loci = loc{k,i};
        if isempty(loci) && sub_plot == 0 % if not making subplots, and there is nothing to plot ...
            close(fig_close); % close the figure window
        end
        for j = 1:length(loci) % run through all scenarios that have data for state 1
            [~,locj] = ismember(loci(j),locall); % to find the right colour ...
            if bw == 1 % plot in black and white
                plot(t,Xall{i}(:,loci(j)),'k-','LineWidth',1) % plot each model curve with a black line
            else
                if locj <= 6
                    plot(t,Xall{i}(:,loci(j)),[plotcol(locj),'-'],'LineWidth',2) % plot each model curve with a different color
                elseif locj <= 12 % if we have more than 6 treatments, re-use colours and make broken line
                    plot(t,Xall{i}(:,loci(j)),[plotcol(locj-6),'--'],'LineWidth',2) % plot each model curve with a different color
                elseif locj <= 18 % if we have more than 12 treatments, re-use colours and make dot-stripe line
                    plot(t,Xall{i}(:,loci(j)),[plotcol(locj-12),'.-'],'LineWidth',2) % plot each model curve with a different color
                else
                    disp('cannot plot more than 18 lines yet ... so modify your script')
                end
            end
            % Note: if sho=0 this only plots scenarios that actually have data!
        end
        
        if bw ~= 1 % plot in color: then use colors for the legend
            % add a legend using the cell array L
            if ~isempty(Le{k,i}) && legsup ~= 1
                h_leg{k,i} = legend(Le{k,i}); % create a legend
                set(h_leg{k,i},ft.name,ft.legend); % use standard font formats
            else
                h_leg{k,i} = [];
            end
        end
        
        plotdata = DATA{k,i};
        plotweights = W{k,i}; % also collect weights ...
        lam = plotdata(1,1);
        % If we set transf to zero, we don't want to modify sub-lethal
        % endpoints with transformations, when calculating means to plot.
        if lam >= 0 && transf == 0
            lam = 1;
        end
        
        if lam == -1 || lam == -3 % if we have survival data ...
            ylim([0 1.02]); % restrict axis from 0-1
            if repls == 1 % if we have survival data and we plot individual replicates
                % Then we need to recalculate survivors into survival
                % probability (taking care of NaNs and missing/removed
                % animals). This is dealt with in the function recalc_data.
                % If we plot means, the transformation is dealt with below,
                % for each set of replicates since they will be combined.
                a = plotdata(2:end,2:end); % extract data set
                for i_recalc = 1:size(a,2) % run through all replicates (recalc_data does one replicate at a time)
                    a_tmp = recalc_data(a(:,i_recalc),plotweights(:,i_recalc),-1); % -1 signals survival data
                    a(:,i_recalc) = a_tmp(:,1); % replace survivors in 'a'  by survival probability (corrected where needed)
                    % survival probability is first column in output of recalc_data
                end
                
                if lam == -3 % we have healthy/immobile/dead states
                    if isfield(glo,'loc_i') % this is needed since we also have a reduced, simplified version!
                        if ismember(i,[glo.loc_i glo.loc_d glo.loc_id]) % than we have an immobility or death state
                            a(1,:) = 0; % make all initial things at t=0 zero!
                            % Needed as first row in data for immobility and death
                            % states will be the total amount of animals with which
                            % the experiment started!
                        end
                    end
                end
                
                plotdata(2:end,2:end) = a; % and put them back
            end
        end
        
        % for survival dose-response curves
        if lam == -2 % if we have survival data in a dose-response setting ...
            a = plotdata(2:end,2:end); % extract data set
            w = W{k,i}; % extract weight factors, which is now the number of start animals
            plotdata(2:end,2:end) = a ./ w; % recalculate to survival prob. and put them back
            ylim([0 1.02]); % restrict axis from 0-1
        end
        
        td = plotdata(2:end,1); % time vector for this data set
        for j = 1:length(loci) % run through all scenarios that have data for state 1
            ind_plot = plotdata(1,:) == X0mat(1,loci(j)); % location to data to plot
            ind_plot(1) = 0; % zero makes sure that time column is not plotted itself

            if sum(ind_plot) > 0 % skip code below in case we simulate and don't fit
                plotpts = plotdata(2:end,ind_plot); % data vector or matrix (if there are replicates) for treatment j
                plotwts = plotweights(:,ind_plot(2:end));
                
                if repls == 2 % make boxplots rather than plotting just means or replicates
                    boxplot(plotpts',td,'positions',td,'PlotStyle','compact');%,'Colors','g'); % ,'Notch','on')
                    % boxplot(plotpts',td,'positions',td,'PlotStyle','traditional','Colors','k','Notch','on')
                    
                    % boxplot messes up the axes, so return everything to auto!
                    ylim auto
                    xlim auto
                    yticks('auto')
                    xticks('auto')
                    yticklabels('auto')
                    xticklabels('auto')
                    max_data = max(max_data,max(plotpts(:))); % remember maximum of plotted data
                    min_data = min(min_data,min(plotpts(:))); % remember minimum of plotted data
                    plotpts = nan(size(td)); % make plotpts NaN so no regular symbols are plotted
                    plotwts = plotpts; % also make the weights NaN
                end
                
                if repls == 0 % plot means instead of replicates
                    % Use recalc_data to construct means from the
                    % replicates. This is done for survival and sub-lethal
                    % endpoints. For survival, a translation to
                    % probabilities is made. If replicates are plotted,
                    % this translation to probabilities has already been
                    % made for the entire data set above.
                    
                    if lam == -1 % survival data
                        chck_nans = mean(isnan(plotpts),2); % this is 0 when there are no NaNs, and 1 when there are NaNs in all replicates
                        if ~all(chck_nans == 1 | chck_nans == 0) % if they are NOT all 0 or all 1, set a flag as we cannot plot some means!
                            flag_meanskip = 1; % set a flag so we generate a warning on screen
                        end
                    end
                    
                    [a_tmp,~,w_out] = recalc_data(plotpts,plotwts,lam);
                    plotpts = a_tmp(:,1); % replace plot points by combined value (corrected where needed)
                    plotwts = w_out; % sum of weights, with NaNs in data as zero weight
                    % Note: plotwts is only used for plotting outliers on
                    % sub-lethal data below.
                    
                    if lam == -3 % we have healthy/immobile/dead states
                        if ismember(i,[glo.loc_i glo.loc_d glo.loc_id]) % than we have an immobility or death state
                            plotpts(1,:) = 0; % make all initial things at t=0 zero!
                            % Needed as first row in data for immobility and death
                            % states will be the total amount of animals with which
                            % the experiment started!
                        end
                    end
                    
                    % % plot CIs on the means! this rapidly becomes messy
                    % for i_ci = 1:length(td)
                    %     plot([td(i_ci) td(i_ci)],[a_tmp(i_ci,2) a_tmp(i_ci,3)],'k-')
                    % end
                end
                
                max_data = max(max_data,max(plotpts(:))); % remember maximum of plotted data
                min_data = min(min_data,min(plotpts(:))); % remember minimum of plotted data
                
                [~,locj] = ismember(loci(j),locall); % to find the right colour ...
                if bw == 1 % plot in black and white
                    if cn == 1 % plot connection between line and point
                        plotcon = interp1(t,Xall{i}(:,loci(j)),td); % model prediction at measured data points
                        for con = 1:length(td)
                            for con2 = 1:size(plotpts,2)
                                if ~isnan(plotpts(con,con2)) % only do that when there is a data point
                                    % hc{j}(con,con2) = plot([td(con) td(con)],[plotcon(con) plotpts(con,con2)],'k:','LineWidth',1);
                                    plot([td(con) td(con)],[plotcon(con) plotpts(con,con2)],'k:','LineWidth',1);
                                end
                            end
                        end
                    end
                    
                    if locj <= 8 % for BW plots, handle h is used for the legend, h_top for putting symbols on top
                        htmp = plot(td,plotpts,['k',markers_dr{locj}],'MarkerFaceColor','w');
                    elseif locj <= 16
                        htmp = plot(td,plotpts,['k',markers_dr{locj-8}],'MarkerFaceColor','k');
                    elseif locj <= 24
                        htmp = plot(td,plotpts,['k',markers_dr{locj-16}],'MarkerFaceColor','y');
                    end
                    h(j) = htmp(1); % needed for replicates, to get 1 handle number for the legend

                    % collect the handles for all symbols, to later put symbols on top
                    for i_top = 1:length(htmp)
                        h_top_tmp = [];
                        h_top_tmp(1) = htmp(i_top); % this forces the handle to a number!
                        h_top = [h_top h_top_tmp];
                    end
                    
                else % for color plots, handle h_top is used only for putting symbols on top
                    if locj <= 6
                        htmp = plot(td,plotpts,'ko','MarkerFaceColor',plotcol(locj));
                    elseif locj <= 12 % if we have more than 6 treatments, re-use colours and make squares
                        htmp = plot(td,plotpts,'ks','MarkerFaceColor',plotcol(locj-6));
                    elseif locj <= 18 % if we have more than 12 treatments, re-use colours and make triangles
                        htmp = plot(td,plotpts,'k^','MarkerFaceColor',plotcol(locj-12));
                    end
                    for i_top = 1:length(htmp) % there might be more than 1 handle if there are replicates
                        h_top_tmp = [];
                        h_top_tmp(1) = htmp(i_top); % this forces the handle to a number!
                        h_top = [h_top h_top_tmp];
                    end
                end
                
                % TEST identify outliers (points with weight zero) while plotting
                if lam ~= -1 && lam ~= -3 && any(plotwts(:)==0) && outl == 1 % if we do NOT have survival data ...
                    plot_test = plotpts; % copy plot points
                    plot_test(plotwts>0) = NaN; % remove points with weight>zero (so only outlier are left)
                    
                    % plot outliers with a specific symbol
                    if bw == 1 % plot in black and white
                        htmp_tmp = plot(td,plot_test,['k',markers_dr{locj}],'MarkerFaceColor',[0.7 0.7 0.7]);
                    else
                        htmp_tmp = plot(td,plot_test,'ko','MarkerSize',10);
                        % for now, only do this for BW ... for color a
                        % bigger marker ...
                        % htmp_tmp = [];
                    end
                    % for now, plot the weight zero data as black
                    % symbols; this is not handy when there are more
                    % than 8 scenarios ...
                    %
                    % collect the handles for outlier symbols, to later put symbols on top
                    for i_top = 1:length(htmp_tmp)
                        h_out_tmp = [];
                        h_out_tmp(1) = htmp_tmp(i_top); % this forces the handle to a number!
                        h_out = [h_out h_out_tmp];
                    end
                end
                
            end
        end
        
        if limax == 1 && ~isempty(plotpts) && numel(plotpts) > 1
            % axis([min(td) max(td)*1.02 min_data/1.02 max_data*1.02]); % limit both axes to data
            xlim([min(td) max(td)*1.02]); % limit only x-axes to data
        else
            if ~isempty(td)
                max_x = max(max(t),max(td));
                min_x = min(min(t),min(td));
            else
                max_x = max(t);
                min_x = min(t);
            end
            diff_x = 0.01*(max_x-min_x);
            if min_x ~= 0
                xlim([min_x-diff_x max_x+diff_x])
            else
                xlim([0 max_x+diff_x])
            end
        end
        
        if exist('Xlo','var') % do we have max and min model lines from the slice sampling?
 
            for j = 1:length(loci) % run through all scenarios that have data for state i
                
                if bw == 1 % plot in black and white
                    plot(t,Xlo{i}(:,loci(j)),'k:','LineWidth',1) % plot each confidence line as broken line
                    plot(t,Xhi{i}(:,loci(j)),'k:','LineWidth',1) % plot each confidence line as broken line
                else
                    [~,locj] = ismember(loci(j),locall); % to find the right colour ...
                    if locj <= 6
                        plot(t,Xlo{i}(:,loci(j)),[plotcol(locj),':'],'LineWidth',2) % plot each model curve with a different color
                        plot(t,Xhi{i}(:,loci(j)),[plotcol(locj),':'],'LineWidth',2) % plot each model curve with a different color
                    elseif locj <= 12 % if we have more than 6 treatments, re-use colours
                        plot(t,Xlo{i}(:,loci(j)),[plotcol(locj-6),':'],'LineWidth',2) % plot each model curve with a different color
                        plot(t,Xhi{i}(:,loci(j)),[plotcol(locj-6),':'],'LineWidth',2) % plot each model curve with a different color
                    elseif locj <= 18 % if we have more than 12 treatments, re-use colours
                        plot(t,Xlo{i}(:,loci(j)),[plotcol(locj-12),':'],'LineWidth',2) % plot each model curve with a different color
                        plot(t,Xhi{i}(:,loci(j)),[plotcol(locj-12),':'],'LineWidth',2) % plot each model curve with a different color
                    end
                end
            end
            
            % Repeat that for the bounds including sampling error, for survival data only!
            if isfield(glo,'locS') && i == glo.locS && numel(XloS)>1 % do we have max and min model lines from the sampling error calc?
                for j = 1:length(loci) % run through all scenarios that have data for state i
                    
                    if bw == 1 % plot in black and white
                        plot(t,XloS(:,loci(j)),'k--','LineWidth',0.5) % plot each model curve as broken line
                        plot(t,XhiS(:,loci(j)),'k--','LineWidth',0.5) % plot each model curve as broken line
                    else
                        
                        [~,locj] = ismember(loci(j),locall); % to find the right colour ...
                        if locj <= 6
                            plot(t,XloS(:,loci(j)),[plotcol(locj),'--'],'LineWidth',1) % plot each model curve with a different color
                            plot(t,XhiS(:,loci(j)),[plotcol(locj),'--'],'LineWidth',1) % plot each model curve with a different color
                        elseif locj <= 12 % if we have more than 6 treatments, re-use colours
                            plot(t,XloS(:,loci(j)),[plotcol(locj-6),'--'],'LineWidth',1) % plot each model curve with a different color
                            plot(t,XhiS(:,loci(j)),[plotcol(locj-6),'--'],'LineWidth',1) % plot each model curve with a different color
                        elseif locj <= 18 % if we have more than 12 treatments, re-use colours
                            plot(t,XloS(:,loci(j)),[plotcol(locj-12),'--'],'LineWidth',1) % plot each model curve with a different color
                            plot(t,XhiS(:,loci(j)),[plotcol(locj-12),'--'],'LineWidth',1) % plot each model curve with a different color
                        end
                    end
                end
            end
        
            uistack(h_top,'top') % put symbols on top, at the end of all plotting
            uistack(h_out,'top') % put outlier symbols on top, at the end of all plotting
        end
       
        if y_zero == 1 % then we modify the y-axis to always start at zero
            YlimCurrent = get(gca,'YLim');
            ylim([0 YlimCurrent(2)])
        end
        
        if bw == 1 % plot in black and white: then use symbols for the legend
            % add a legend using the cell array L
            if ~isempty(Le{k,i}) && ~isempty(plotpts) && legsup ~= 1 % don't show legend if there's no data plotted
                h_leg{k,i} = legend(h(h>0),Le{k,i}(find(h>0))); % create a legend
                set(h_leg{k,i},ft.name,ft.legend); % use standard font formats
            else
                h_leg{k,i} = [];
            end
        end
        
        % option to plot graph on log-scale for dose-response analysis
        if isfield(glo,'logsc') && glo.logsc == 1
            if td(1) == 0 % if the control is really zero
                if isfield(glo,'logzero') && glo.logzero ~= 0 % override default plotting for zero
                    td0 = glo.logzero;
                else
                    td0 = 10^(log10(td(2)) - 2*mean(diff(log10(td(td>0)))));
                end
                plot(td0,plotpts(1,:),'ko','MarkerFaceColor','w') % plot it at a reasonable place as open symbol
            end
            h_ax = gca;
            h_ax.XScale = 'log'; % put x-axis on log-scale
            xlim([td0/1.2 td(end)*1.2])
        end
        
%         % TEST plot events in there as well, using ylim to determine range on y-axis
%         plot([tt{1} tt{1}],ylim.*[1 0.99],'k:','LineWidth',1.5)
        
        if sub_plot == 0 && glo.saveplt > 0 % if we want to save the plot, and individual plots are made
            if n_D > 1
                savenm = ['fit_',savenm1,'_',num2str(k),'_',num2str(i)]; % name includes state and data set
            else
                savenm = ['fit_',savenm1,'_',num2str(i)]; % name only includes state as there is 1 data set
            end
            save_plot(gcf,savenm);
        end
    end
    
    % Making an additional sub-plot with the parameter estimates
    if annot == 1 && sub_plot == 1 % make an extra subplot with the parameter estimates
        h_annot = subplot(m,n,n_X+1); % make an extra subplot, and remember the handle
        h1 = gca; % remember the current axis
        dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
        str_annot = []; % clear the string if it already exists, and fill it step by step
        
        if ~isfield(par_plot,'tag_fitted')
            pmat_plot(:,2) = 0; % set all parameters to 'not fitted'
            str_annot{1} = 'Parameter estimates (not fitted)';
        else
            str_annot{1} = 'Parameter estimates (fitted)';
        end
        str_annot{2} = '=====================================';
        nfields = length(names); % use the names from the global
        for i = 1:nfields % display all parameters on plot
            if pmat_plot(i,5) == 0
                str_annot{2+i} = sprintf('%-6s %10.4g (fit: %1.0f) log-scale',names{i},pmat_plot(i,1),pmat_plot(i,2));
            else
                str_annot{2+i} = sprintf('%-6s %10.4g (fit: %1.0f)',names{i},pmat_plot(i,1),pmat_plot(i,2));
            end
        end
        str_annot{3+i} = '=====================================';
        str_annot{4+i} = sprintf('Min log-lik.: %1.6g (AIC=%1.6g)',FVAL,AIC);
        str_annot{5+i} = sprintf('Filename: %s',glo.basenm); % name of the mydata file
        ti = clock; % what is the current time?
        str_annot{6+i} = sprintf('Analysis date: %s (%02.0f:%02.0f)',date,ti(4),ti(5));
        
        if isfield(glo,'str_extra')
            str_extra = glo.str_extra;
            str_annot{7+i} = '=====================================';
            for i_ex = 1:length(str_extra)
                str_annot(end+1) = str_extra(i_ex);
            end
            str_annot{end+1} = '=====================================';
        end
        
        % make a textbox and place it in the additional subplot
        annotation('textbox',dim,'String',str_annot,'FitBoxToText','on','BackGroundColor','w','Interpreter','none',ft.name,ft.annot);
        delete(h1) % throw away the axes of the plot, the text box is enough
    end
    
    if annot == 2 && sub_plot == 1 % make an extra subplot with the legend
       % PM this needs to be fine-tuned!
       h_annot = subplot(m,n,n_X_plot+1); % make an extra subplot, and remember the handle
       h1 = gca; % remember the current axis
       set(h1,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
       dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
       hold on
       
%        for j = 1:length(L) % run through complete set of legend entries
%            if bw == 1
%                if j <= 8 % for BW plots, handle h is used for the legend, h_top for putting symbols on top
%                    plot(0,0,['k',markers_dr{j}],'MarkerFaceColor','w');
%                elseif j <= 16
%                    plot(0,0,['k',markers_dr{j-8}],'MarkerFaceColor','k');
%                elseif j <= 24
%                    plot(0,0,['k',markers_dr{j-16}],'MarkerFaceColor','y');
%                end 
%            else
%                if j <= 6
%                    plot(0,0,[plotcol(j),'-'],'LineWidth',2) % plot each model curve with a different color
%                elseif j <= 12 % if we have more than 6 treatments, re-use colours and make broken line
%                    plot(0,0,[plotcol(j-6),'--'],'LineWidth',2) % plot each model curve with a different color
%                elseif j <= 18 % if we have more than 12 treatments, re-use colours and make dot-stripe line
%                    plot(0,0,[plotcol(j-12),'.-'],'LineWidth',2) % plot each model curve with a different color
%                end
%            end
%        end
       
       % first find all plotted lines in this plot window (this is a
       % bit complex as there may be supressed states/scenarios!)
       loc_all = [];
       for j = 1:length(loc(k,:)) % run through states that are plotted
           loc_all = [loc_all loc{k,j}];
       end
       loc_all = unique(loc_all);

       % Modified to limit the legend to the scenarios that are shown 
       for j = 1:length(loc_all) % run through complete set of legend entries
           loc_j = loc_all(j);
           if bw == 1
               if loc_j <= 8 % for BW plots, handle h is used for the legend, h_top for putting symbols on top
                   plot(0,0,['k',markers_dr{loc_j}],'MarkerFaceColor','w');
               elseif loc_j <= 16
                   plot(0,0,['k',markers_dr{loc_j-8}],'MarkerFaceColor','k');
               elseif loc_j <= 24
                   plot(0,0,['k',markers_dr{loc_j-16}],'MarkerFaceColor','y');
               end 
           else
               if loc_j <= 6
                   plot(0,0,[plotcol(loc_j),'-'],'LineWidth',2) % plot each model curve with a different color
               elseif loc_j <= 12 % if we have more than 6 treatments, re-use colours and make broken line
                   plot(0,0,[plotcol(loc_j-6),'--'],'LineWidth',2) % plot each model curve with a different color
               elseif loc_j <= 18 % if we have more than 12 treatments, re-use colours and make dot-stripe line
                   plot(0,0,[plotcol(loc_j-12),'.-'],'LineWidth',2) % plot each model curve with a different color
               end
           end
       end
       
       h_leg_tot = legend(L(loc_all)); % create a legend with selected entries (use state 1 for that)
       set(h_leg_tot,ft.name,ft.legend); % use standard font formats

       xlim([1 10]) % make sure the 'data' are off screen, so invisible
       
       dimleg = h_leg_tot.Position; % position of LAST legend made
       % move it to the new subplot and set it to top-left of panel
       h_leg_tot.Position = [dim(1) dim(2)+dim(4)-dimleg(4) dimleg(3:4)];

       h1.Position = [dim(1) dim(2)+dim(4)-dimleg(4) 0.01 0.01];
       
       for j_tmp = 1:n_X_plot % delete the other legends
           j = show_X(j_tmp);
           h_leg_tmp = h_leg{k,j}; % take legend handle
           delete(h_leg_tmp); % and delete it
       end 
       h1.Visible = 'off'; % hide axes
    end
    
    if sub_plot == 1
        if notitle == 0
            if ~exist('Xlo','var') % if making sub-plots, put a single title on top of the graph
                % if there are CIs, you get a different title
                ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
                if n_D == 1
                    h_txt{k} = text(0.5, 1,['Plotted from: ',savenm1, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
                else
                    h_txt{k} = text(0.5, 1,['Plotted from: ',savenm1, ' (',date,'), series ',num2str(k),' of ',num2str(n_D)],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
                end
                set(h_txt{k},ft.name,ft.text); % use standard formatting for this header
            end
            if exist('Xlo','var') % if making sub-plots, put a single title on top of the graph
                % if there are CIs, you get a different title
                ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
                switch type_conf
                    case 1
                        h_txt{k} = text(0.5, 1,'Intervals: 95% of model predictions from posterior','HorizontalAlignment','center','VerticalAlignment', 'top');
                    case 2
                        h_txt{k} = text(0.5, 1,'Intervals: predictions from likelihood region (likreg)','HorizontalAlignment','center','VerticalAlignment', 'top');
                    case 3
                        h_txt{k} = text(0.5, 1,'Intervals: predictions from likelihood region (parspace)','HorizontalAlignment','center','VerticalAlignment', 'top');
                end
                set(h_txt{k},ft.name,ft.text); % use standard formatting for this header
            end
        end
        
        if glo.saveplt > 0 % if we want to save the plot
            savenm = ['fit_',savenm1,'_',num2str(k)];%
            save_plot(gcf,savenm,h_txt{k});
        end
    end
       
end

if flag_meanskip == 1 && repls == 0 % we could not calculate some means, and plot means
    warning('off','backtrace') % no need to display where the warning is generated
    warning('For survival, some means are not plotted because there are NaNs in some replicates.')
    warning('on','backtrace'), disp(' ')
end

%% Plot zero-variate data points

if isfield(par_plot,'tag_fitted') && ~exist('Xlo','var') % only when fitted and not called by calc_conf
    if zvd_plot == 1 % Plot zero-variate data in a nice way ...
        if ~isempty(zvd) && ~isempty(namesz)
            make_fig(1,1); % create figure window of correct size
            for i = 1:length(namesz) % go through contents of the global zvd
                a = eval(['zvd.',namesz{i}]); % extract zero-variate info from zvd_out
                errorbar(i-0.1,a(1),1.96*a(2),1.96*a(2),'ko','MarkerFaceColor','w')
                hold on
                plot(i+0.1,a(3),'ko','MarkerFaceColor','y')
                if i == 1
                    h_leg = legend('zero-var. data with CI','model result');
                    set(h_leg,ft.name,ft.legend); % use standard font formats
                end
            end
            xlim([0 i+1]);
            ym = ylim; % ask for current y limits
            ylim([0 ym(2)*1.02]); % start at zero, but keep upper limit
            
            set(gca,'xtick',1:i,'xticklabel',namesz,'LineWidth',1,ft.name,ft.ticks)
            xlabel('zero-variate data point name',ft.name,ft.label)
            ylabel('zero-variate value',ft.name,ft.label)
            
            if glo.saveplt > 0 % if we want to save the plot
                savenm = ['zerovar_',savenm];%
                save_plot(gcf,savenm);
            end
        end
    end
end

%% Plot extra data points on different x-axes
% Note: the Calanus DEBkiss model makes use of this option

if n_X2 > 0 % if there are additional data fitted ...
    plotcol = 'krbgcm'; % 6 colors should be enough ...
    markers_dr = {'s','^','o','d','v','p','>','<'}; % complete set of 8 markers

    % make one plot with subplots for each data set
    n = ceil(sqrt(n_X2));
    m = ceil(n_X2/n);
    
    if sub_plot == 1
        make_fig(m,n); % create a figure window of correct size
    end

    for i = 1:n_X2 % run through extra data sets
        h = [];
        
        % make a subplot for set i
        if sub_plot == 1 % make a subplot for state i
            subplot(m,n,i)
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        else
            make_fig(1,1);  % create a figure window of correct size
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        end
        hold on
        
        % add labels to the axes
        xlabel(glo.xlab2{i},ft.name,ft.label)
        ylabel(glo.ylab2{i},ft.name,ft.label)
        
        if notitle == 0
            % add a title to the graph with the filename in it
            title(['Plotted from: ',savenm1, ' (',date,')'],'Interpreter','none',ft.name,ft.title)
        end
          
        [~,loci]=ismember(DATAx{i}(1,2:end),X0mat(1,:)); % loci is location of state i in the scenario range
        % loci = X0mat(1,:);
        
        loci(loci==0) = []; % remove zeros (where data for certain scenarios are ignored)
        loci = unique(loci); % only keep the unique entries (in case of replicated data)
        
        for j = 1:length(loci) % run through all scenarios that have data for state 1
            [~,locj] = ismember(loci(j),locall); % to find the right colour ...
            if bw == 1 % plot in black and white
                plot(Xall2x{i}(:,loci(j)),Xall2y{i}(:,loci(j)),'k-','LineWidth',1) % plot each model curve with a black line
            else
                if locj <= 6
                    h(j) = plot(Xall2x{i}(:,loci(j)),Xall2y{i}(:,loci(j)),[plotcol(locj),'-'],'LineWidth',2); % plot each model curve with a different color
                elseif locj <= 12 % if we have more than 6 treatments, re-use colours and make broken line
                    h(j) = plot(Xall2x{i}(:,loci(j)),Xall2y{i}(:,loci(j)),[plotcol(locj-6),'--'],'LineWidth',2); % plot each model curve with a different color
                elseif locj <= 18 % if we have more than 12 treatments, re-use colours and make dot-stripe line
                    h(j) = plot(Xall2x{i}(:,loci(j)),Xall2y{i}(:,loci(j)),[plotcol(locj-12),'.-'],'LineWidth',2); % plot each model curve with a different color
                else
                    disp('cannot plot more than 18 lines yet ... so modify your script')
                end
            end
        end
                
        for j = 1:length(loci) % run through all scenarios that have data for state 1
            
            % if ismember(loci(j),DATAx{i}(1,2:end))
            
            [~,locj] = ismember(loci(j),locall); % to find the right colour ...
            locx = find(X0mat(1,loci(j)) == DATAx{i}(1,2:end)); % this should work if we have more columns of data in one scenario
            
            plotdata = DATAx{i}(2:end,[1 locx+1]);
            plotoutl = plotdata; % make a copy for the outliers
            plotdata(Wx{i}(:,locx) == 0,2) = NaN; % make outlier NaN
            plotoutl(Wx{i}(:,locx) > 0,2) = NaN; % make non-outlier NaN
                        
            if bw == 1 % plot in black and white
                
                if cn == 1 % plot connection between line and point
                    plotcon = interp1(Xall2x{i}(:,loci(j)),Xall2y{i}(:,loci(j)),DATAx{i}(2:end,1)); % model prediction at measured data points
                    for con = 1:length(DATAx{i}(2:end,1))
                        if ~isnan(DATAx{i}(con+1,locx+1)) % only do that when there is a data point
                            plot([DATAx{i}(con+1,1) DATAx{i}(con+1,1)],[plotcon(con) DATAx{i}(con+1,locx+1)],'k:','LineWidth',1);
                        end
                    end
                end
                
                % DATAx{i}(2:end,1),DATAx{i}(2:end,locx+1)
                if locj <= 8
                    h_tmp = plot(plotdata(:,1),plotdata(:,2),['k',markers_dr{locj}],'MarkerFaceColor','w');
                elseif locj <= 16
                    h_tmp = plot(plotdata(:,1),plotdata(:,2),['k',markers_dr{locj-8}],'MarkerFaceColor','k');
                elseif locj <= 24
                    h_tmp = plot(plotdata(:,1),plotdata(:,2),['k',markers_dr{locj-16}],'MarkerFaceColor','y');
                end
                plot(plotoutl(:,1),plotoutl(:,2),['k',markers_dr{locj}],'MarkerFaceColor',[0.7 0.7 0.7]);
                h(j) = h_tmp(1); % this is needed in case there are multiple columns of extra data for one scenario
                
            else
                if locj <= 6
                    plot(plotdata(:,1),plotdata(:,2),'ko','MarkerFaceColor',plotcol(locj));
                    plot(plotoutl(:,1),plotoutl(:,2),'ko','MarkerFaceColor',plotcol(locj));
                    plot(plotoutl(:,1),plotoutl(:,2),'ko','MarkerSize',10);
                elseif locj <= 12 % if we have more than 6 treatments, re-use colours and make squares
                    plot(plotdata(:,1),plotdata(:,2),'ks','MarkerFaceColor',plotcol(locj-6));
                    plot(plotoutl(:,1),plotoutl(:,2),'ko','MarkerFaceColor',plotcol(locj-6));
                    plot(plotoutl(:,1),plotoutl(:,2),'ko','MarkerSize',10);
                elseif locj <= 18 % if we have more than 12 treatments, re-use colours and make triangles
                    plot(plotdata(:,1),plotdata(:,2),'k^','MarkerFaceColor',plotcol(locj-12));
                    plot(plotoutl(:,1),plotoutl(:,2),'ko','MarkerFaceColor',plotcol(locj-12));
                    plot(plotoutl(:,1),plotoutl(:,2),'ko','MarkerSize',10);
                end
                
            end
            % end
        end
        
        % add a legend using the cell array L
        if ~isempty(h) && legsup ~= 1
            h_leg = legend(h,L{loci}); % create a legend TEST
            set(h_leg,ft.name,ft.legend); % use standard font formats
        end
        
    end
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['extra_',savenm];%
        save_plot(gcf,savenm);
    end
    
end

glo.h_txt = h_txt; % make the text handles available in global, for suppressing them in the script, if needed
