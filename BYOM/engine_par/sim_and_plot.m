function figh = sim_and_plot(par,opt_sim_input)

% Function to simulate and plot (systems of) ODEs 'real time'
% Various plotting options are available
% 1) 3d-plot (specify which 3 states, in which order, in your byom script)
% 2) 2d-plot (specify which 2 states, in which order, in your byom script)
% 3) subplots for each state (scenarios as different color)
% 4) subplots for each scenario (states are different color)
% 5) dx/dt versus x for each state in a separate subplot
% 6) 2d plot for each scenario separate
% 7) 3d plot for each scenario separate
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo X0mat

ylab = glo.ylab;
xlab = glo.xlab;
leglab1 = glo.leglab1;
leglab2 = glo.leglab2;
t = glo.t;

c  = X0mat(1,:);      % scenarios
ns = size(X0mat,1)-1; % number of states
nc = length(c);       % number of scenarios

plotcol  = 'krbgcmy'; % 7 colors should be enough ...
trainsz  = 8; % plot size of the running 'train'
trailsz  = 1; % plot size of the permanent 'trail' (zero for NO trail)
plotstrt = 1; % plot starting point in 2d and 3d plots
plotvct  = 0; % plot field of vector arrows in 2d plots when possible

% =========================================================================
% Make sure that versions of Matlab 2017a and newer do not automatically
% update figure legends when new info is plotted in the same plot. This
% cannot be set for older versions (they don't update the legend anyway).
% I first modified the legend statements in the plotting routines to
% prevent updating, but this yields errors with older versions.
if ~verLessThan('matlab','9.2')
    % -- Code to run in MATLAB R2017a and later here --
    set(groot,'defaultLegendAutoUpdate','off');
end
% This is also done in prelim_checks, but the simulation scripts will
% usually skip that script entirely!
% =========================================================================

% Options for simulations with sim_and_plot, with default settings. Need to
% do this differently for sim_and_plot as prelim_checks will be skipped!
opt_sim.plottype = 4;       % type of simulation plot to make (see sim_and_plot)
opt_sim.stxyz    = [1 2 3]; % states to plot on x, y and z-axis
opt_sim.plot_int = 10;      % interval for plotting (how many time points to do in one go)
opt_sim.axrange  = [];      % axis ranges for each state, default: autoscale
opt_sim.Xeq      = [];      % equilibrium points to plot
opt_sim.ylog     = 0;       % set to 1 to put y-axis on log scale (for type 3-4)
opt_sim.yzero    = 0;       % set to 1 to minimise y-axis to zero (for type 3-4)

names=fieldnames(opt_sim);
for i = 1:length(names)
    if isfield(opt_sim_input,names{i}) % option is in the input structure
        eval([names{i},'=opt_sim_input.',names{i},';'])
    else % use the default
        eval([names{i},'=opt_sim.',names{i},';'])
    end
end

stx = stxyz(1);
sty = stxyz(2);
stz = stxyz(3);

plot_int = max(3,plot_int); % less than 3 gives problems with calling call_deri

% some preliminary checks ...
if ns > 7 && plottype == 4
    error('Too many states for this plottype')
end
if nc > 7 && plottype == 3
    error('Too many scenarios for this plottype')
end
if ns < 3 && plottype == 1
    error('Not enough state variable for a 3d plot')
end
if ns < 2 && (plottype == 2 || plottype == 6)
    error('Not enough state variable for a 2d plot')
end
if ~exist('ylab','var') || isempty(ylab{1})
    ylab{1} = 'State var. 1';
end
if ns ~= size(ylab,2)
    for i = 2:ns
        if length(ylab) < i
            ylab{i} = ['State var. ',num2str(i)];
        end
    end
end
if ~exist('xlab','var')
    xlab = 'time';
end
if ~exist('leglab1','var')
    leglab1 = 'Scen. ';
end
if ~exist('leglab2','var')
    leglab2 = '';
end
% =========================================================================

switch plottype
    case 3
        n = ceil(sqrt(ns));
        m = ceil(ns/n);
        [figh,ft] = make_fig(m,n);
    case 4
        n = ceil(sqrt(nc));
        m = ceil(nc/n);
        [figh,ft] = make_fig(m,n);
    case 5
        n = ceil(sqrt(ns));
        m = ceil(ns/n);
        [figh,ft] = make_fig(m,n);
    case 6
        n = ceil(sqrt(nc));
        m = ceil(nc/n);
        [figh,ft] = make_fig(m,n);
    case 7
        n = ceil(sqrt(nc));
        m = ceil(nc/n);
        [figh,ft] = make_fig(m,n);
    otherwise
        [figh,ft] = make_fig(2,2); % there is only one plot, but a big one is nice!
end

hold on
axis square
figure(figh); % bring figure to front

% % TESTING VIDEO MAKING
% v = VideoWriter('test_sim.avi');
% v.FrameRate = 20; %  frames per second
% open(v);

switch plottype
    case 1
        
        % =================================================================
        % 3D plot
        % =================================================================
        view(-20,25) % camera angle
        grid on
        xlabel(ylab{stx},ft.name,ft.label)
        ylabel(ylab{sty},ft.name,ft.label)
        zlabel(ylab{stz},ft.name,ft.label)
        title(['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter', 'none',ft.name,ft.title)
        set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        
        if ~isempty(axrange) % if there is an axis range, use it
            axis([axrange(stx,:) axrange(sty,:) axrange(stz,:)])
        end
        if ~isempty(Xeq) % if equilibria are given, plot them
            for ieq = 1:size(Xeq,1)
                stem3(Xeq(ieq,stx),Xeq(ieq,sty),Xeq(ieq,stz),'ko','MarkerFaceColor','y','BaseValue',axrange(stz,1),'LineWidth',2)
                % plot3(Xeq(ieq,stx),Xeq(ieq,sty),Xeq(ieq,stz),'ko','MarkerFaceColor','y')
            end
        end
        
        X0tmp = X0mat;
        Tn = t(1:plot_int);
        
        % First run to get started
        L = cell(nc,1); % pre-define cell array
        for k = 1:nc % run through scenarios
            if plotstrt == 1
                X0 = X0mat(2:end,k); % extract starting values for this scenario
                plot3(X0(stx),X0(sty),X0(stz),[plotcol(k),'s'],'MarkerFaceColor',plotcol(k)) % plot starting point
            end
            X = call_deri(Tn,par,X0mat(:,k),glo); % use call_deri.m to provide the output
            hl1(k)=plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trainsz);
            if trailsz > 0
                plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trailsz);
            end
            drawnow
            X0tmp(2:end,k) = (X(end,:))';
            L{k} = [leglab1,num2str(X0mat(1,k)),leglab2]; % create legend entries
        end
        h_leg = legend(hl1,L); % plot legend using plot handles
        set(h_leg,ft.name,ft.legend); % use standard font formats
        pause(0.05) % strangely enough, otherwise the legend is ruined ...
        
        % Continue from the first run for the remainder of the time vector
        for i = 2:floor(length(t)/plot_int) % note that some points might not be plotted now!
            for k = 1:nc
                Tn = t((i-1)*plot_int:i*plot_int); % new time vector
                X = call_deri(Tn,par,X0tmp(:,k),glo); % use call_deri.m to provide the output
                delete(hl1(k)) % delete old plot segment before plotting a new one
                hl1(k) = plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trainsz);
                if trailsz > 0
                    plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trailsz);
                end
                drawnow
                X0tmp(2:end,k) = (X(end,:))'; % use previous last as new first
                
                % % TESTING VIDEO MAKING
                % frame = getframe(gcf);
                % writeVideo(v,frame);
                
            end
        end
        if trailsz > 0
            for k = 1:nc
                delete(hl1(k)) % delete old plot segment at the end
            end
        end
        
        % % TESTING VIDEO MAKING
        % close(v);
        
    case 2
        
        % =================================================================
        % 2D plot
        % =================================================================
        
        if ~isempty(Xeq)
            for ieq = 1:size(Xeq,1)
                plot(Xeq(ieq,stx),Xeq(ieq,sty),'ko','MarkerFaceColor','y')
            end
        end
        xlabel(ylab(stx),ft.name,ft.label)
        ylabel(ylab{sty},ft.name,ft.label)
        title(['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter', 'none',ft.name,ft.title)
        set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        
        X0tmp = X0mat;
        Tn = t(1:plot_int);
        
        if ~isempty(axrange) % if there is an axis range, use it
            axis([axrange(stx,:) axrange(sty,:)])
            
            % if there is an axis range, we can also make a vector field!
            if plotvct == 1 && ns == 2 && nc == 1 % only when there are 2 states and 1 scenario
                % if there are more scenarios in one plot, this is too messy
                % (I cannot plot more than 1 vector field in a plot)
                x_grid = linspace(axrange(stx,1),axrange(stx,2),15);
                y_grid = linspace(axrange(sty,1),axrange(sty,2),15);
                dX = nan(length(x_grid),length(y_grid)); % pre-define with nans
                dY = nan(length(x_grid),length(y_grid)); % pre-define with nans
                for i =1:length(x_grid)
                    for j = 1:length(y_grid)
                        X0tmp(stx) = x_grid(i);
                        X0tmp(sty) = y_grid(j);
                        dXtmp = derivatives(0,X0tmp,par,c(1)); % use derivatives to calculate dX/dt
                        dX(i,j) = dXtmp(stx);
                        dY(i,j) = dXtmp(sty);
                    end
                end
                [Xpl,Ypl] = meshgrid(x_grid,y_grid);
                quiver(Xpl,Ypl,dX',dY',1.5)
            end
        end
        
        % First run to get started
        for k = 1:nc
            if plotstrt == 1
                X0 = X0mat(2:end,k);
                plot(X0(stx),X0(sty),[plotcol(k),'s'],'MarkerFaceColor',plotcol(k)) % plot start point
            end
            X = call_deri(Tn,par,X0mat(:,k),glo); % use call_deri.m to provide the output
            hl1(k)=plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trainsz);
            if trailsz > 0
                plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trailsz);
            end
            drawnow
            X0tmp(2:end,k) = (X(end,:))';
            L{k} = [leglab1,num2str(X0mat(1,k)),leglab2]; % create legend entries
        end
        h_leg = legend(hl1,L); % plot legend using plot handles
        set(h_leg,ft.name,ft.legend); % use standard font formats
        
        % Continue from the first run for the remainder of the time vector
        for i = 2:floor(length(t)/plot_int) % note that some points might not be plotted now!
            for k = 1:nc
                Tn = t((i-1)*plot_int:i*plot_int); % new time vector
                X = call_deri(Tn,par,X0tmp(:,k),glo); % use call_deri.m to provide the output
                delete(hl1(k)) % delete old plot segment before plotting a new one
                hl1(k) = plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trainsz);
                if trailsz > 0
                    plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trailsz);
                end
                drawnow
                X0tmp(2:end,k) = (X(end,:))'; % use previous last as new first
            end
        end
        if trailsz > 0
            for k = 1:nc
                delete(hl1(k)) % delete old plot segment at the end
            end
        end
        
    case 3
        
        % =================================================================
        % make one subplot for each state variable
        % =================================================================
        for j = 1:ns % run through states
            h_s{j} = subplot(m,n,j);
            hold on
            xlabel(xlab,ft.name,ft.label)
            ylabel(ylab{j},ft.name,ft.label)
            % title(['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter', 'none',ft.name,ft.title)
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            if ylog == 1
                set(gca,'YScale','log')
            end
            if ~isempty(axrange) % if there is an axis range, use it
                axis([0 t(end) axrange(j,:)])
            else
                xlim([0 t(end)])
            end
        end
        
        X0tmp = X0mat;
        Tn = t(1:plot_int);
        
        ymax = zeros(1,ns); % initialise vector to catch max y
        
        % First run to get started
        L = cell(nc,1); % pre-define cell array
        for k = 1:nc
            X = call_deri(Tn,par,X0mat(:,k),glo); % use call_deri.m to provide the output
            for j = 1:ns
                subplot(h_s{j})
                hold on
                plot(Tn,X(:,j),[plotcol(k),'-'],'LineWidth',2);
                L{k} = [leglab1,num2str(X0mat(1,k)),leglab2];
            end
            
            if ylog == 0 && yzero == 1 % modify y-axis on the fly!
                ymax = max(ymax,max(X,[],1));
                for j = 1:ns
                    ylim(h_s{j},[0 ymax(j)])
                end
            end
            
            drawnow
            X0tmp(2:end,k) = (X(end,:))';
        end
        for j = 1:ns % create legends in each subplot
            subplot(h_s{j})
            h_leg = legend(L,'FontSize',12);
            set(h_leg,ft.name,ft.legend); % use standard font formats
        end
        
        % Continue from the first run for the remainder of the time vector
        for i = 2:floor(length(t)/plot_int) % note that some points might not be plotted now!
            for k = 1:nc
                Tn = t((i-1)*plot_int:i*plot_int); % new time vector
                X = call_deri(Tn,par,X0tmp(:,k),glo); % use call_deri.m to provide the output
                for j = 1:ns
                    subplot(h_s{j})
                    hold on
                    plot(Tn,X(:,j),[plotcol(k),'-'],'LineWidth',2);
                end
                
                if ylog == 0 && yzero == 1 % modify y-axis on the fly!
                    ymax = max(ymax,max(X,[],1));
                    for j = 1:ns
                        ylim(h_s{j},[0 ymax(j)])
                    end
                end
                
                drawnow
                X0tmp(2:end,k) = (X(end,:))'; % use previous last as new first
            end
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
        
    case 4
        
        % =================================================================
        % make one plot for each scenario
        % =================================================================
        for k = 1:nc % run through scenarios
            h_s{k} = subplot(m,n,k);
            hold on
            xlabel(xlab,ft.name,ft.label)
            ylabel(['state variables, ',leglab1,num2str(c(k)),leglab2],ft.name,ft.label)
            % title(['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter', 'none',ft.name,ft.title)
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            if ylog == 1
                set(gca,'YScale','log')
            end
            if ~isempty(axrange) % if there is an axis range, use it
                axis([0 t(end) min(axrange(:)) max(axrange(:))])
            else
                xlim([0 t(end)])
            end
        end
        
        X0tmp = X0mat;
        Tn = t(1:plot_int);
        
        ymax = zeros(1,ns); % initialise vector to catch max y
        
        % First run to get started
        L = cell(ns,1); % pre-define cell array
        for k = 1:nc
            X = call_deri(Tn,par,X0mat(:,k),glo); % use call_deri.m to provide the output
            subplot(h_s{k})
            hold on
            for j = 1:ns
                plot(Tn,X(:,j),[plotcol(j),'-'],'LineWidth',2);
                L{j} = [ylab{j}];
            end
            
            if ylog == 0 && yzero == 1 % modify y-axis on the fly!
                ymax = max(ymax,max(X,[],1));
                for j = 1:ns
                    ylim(h_s{j},[0 ymax(j)])
                end
            end
            
            drawnow
            X0tmp(2:end,k) = (X(end,:))';
        end
        % for k = 1:nc % create legends for each subplot
        %     subplot(m,n,k)
        %     legend(L,'FontSize',12)
        % end
        subplot(m,n,2) % only legend in second plot!
        h_leg = legend(L,'FontSize',12);
        set(h_leg,ft.name,ft.legend); % use standard font formats
        
        % Continue from the first run for the remainder of the time vector
        for i = 2:floor(length(t)/plot_int) % note that some points might not be plotted now!
            for k = 1:nc
                Tn = t((i-1)*plot_int:i*plot_int); % new time vector
                X = call_deri(Tn,par,X0tmp(:,k),glo); % use call_deri.m to provide the output
                subplot(h_s{k})
                hold on
                for j = 1:ns
                    plot(Tn,X(:,j),[plotcol(j),'-'],'LineWidth',2);
                end
                
                if ylog == 0 && yzero == 1 % modify y-axis on the fly!
                    ymax = max(ymax,max(X,[],1));
                    for j = 1:ns
                        ylim(h_s{j},[0 ymax(j)])
                    end
                end
                
                drawnow
                X0tmp(2:end,k) = (X(end,:))'; % use previous last as new first
            end
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
        
    case 5
        
        % =================================================================
        % dx/dt vs x; one subplot for each state variable
        % =================================================================
        for j = 1:ns % run through states
            h_s{j} = subplot(m,n,j);
            hold on
            xlabel(ylab{j},ft.name,ft.label)
            ylabel(['d',ylab{j},'/dt'],ft.name,ft.label)
            % title(['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter', 'none',ft.name,ft.title)
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            
            if ~isempty(axrange) % if there is an axis range, use it
                xlim(axrange(j,:)) % only limit X range
                plot(axrange(j,:),[0 0],'k:')
            end
            if ~isempty(Xeq)
                for ieq = 1:size(Xeq,1)
                    plot(Xeq(ieq,j),0,'ko','MarkerFaceColor','y')
                end
            end
        end
        
        X0tmp = X0mat;
        Tn = t(1:plot_int);
        
        % First run to get started
        L = cell(nc,1); % pre-define cell array
        for k = 1:nc
            X = call_deri(Tn,par,X0mat(:,k),glo); % use call_deri.m to provide the output
            for ki =1:plot_int % run through time vector so that derivatives can stay the same
                dX(ki,:) = (derivatives(Tn(ki),(X(ki,:))',par,c(k)))'; % use derivatives to calculate dX/dt
            end
            for j = 1:ns
                subplot(h_s{j})
                hold on
                hl1(k,j) = plot(X(:,j),dX(:,j),[plotcol(k),'-'],'LineWidth',2);
                L{k} = [leglab1,num2str(X0mat(1,k)),leglab2];
            end
            drawnow
            X0tmp(2:end,k) = (X(end,:))';
        end
        for j = 1:ns % create legends in each subplot
            subplot(h_s{j})
            h_leg = legend(hl1(:,j),L,'FontSize',12);
            set(h_leg,ft.name,ft.legend); % use standard font formats
        end
        
        % Continue from the first run for the remainder of the time vector
        for i = 2:floor(length(t)/plot_int) % note that some points might not be plotted now!
            for k = 1:nc
                Tn = t((i-1)*plot_int:i*plot_int); % new time vector
                X = call_deri(Tn,par,X0tmp(:,k),glo); % use call_deri.m to provide the output
                for ki =1:plot_int+1 % run through time vector so that derivatives can stay the same
                    dX(ki,:) = (derivatives(Tn(ki),(X(ki,:))',par,c(k)))'; % use derivatives to calculate dX/dt
                end
                
                for j = 1:ns
                    subplot(h_s{j})
                    hold on
                    plot(X(:,j),dX(:,j),[plotcol(k),'-'],'LineWidth',2);
                end
                drawnow
                X0tmp(2:end,k) = (X(end,:))'; % use previous last as new first
            end
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
        
    case 6
        
        % =================================================================
        % make one 2d plot for each scenario
        % =================================================================
        for k = 1:nc % run through scenarios
            h_s{k} = subplot(m,n,k);
            hold on
            xlabel(ylab(stx),ft.name,ft.label)
            ylabel(ylab{sty},ft.name,ft.label)
            % title(['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter', 'none',ft.name,ft.title)
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting

            if ~isempty(Xeq)
                for ieq = 1:size(Xeq,1)
                    plot(Xeq(ieq,stx),Xeq(ieq,sty),'ko','MarkerFaceColor','y')
                end
            end
            
            if ~isempty(axrange) % if there is an axis range, use it
                axis([axrange(stx,:) axrange(sty,:)])
                % if there is an axis range, we can also make a vector field!
                if plotvct == 1 && ns == 2 % only when there are 2 states
                    x_grid = linspace(axrange(stx,1),axrange(stx,2),15);
                    y_grid = linspace(axrange(sty,1),axrange(sty,2),15);
                    dX = nan(length(x_grid),length(y_grid)); % pre-define with nans
                    dY = nan(length(x_grid),length(y_grid)); % pre-define with nans
                    for i =1:length(x_grid)
                        for j = 1:length(y_grid)
                            X0tmp(stx) = x_grid(i);
                            X0tmp(sty) = y_grid(j);
                            dXtmp = derivatives(0,X0tmp,par,c(k)); % use derivatives to calculate dX/dt
                            dX(i,j) = dXtmp(stx);
                            dY(i,j) = dXtmp(sty);
                        end
                    end
                    [Xpl,Ypl] = meshgrid(x_grid,y_grid);
                    quiver(Xpl,Ypl,dX',dY',1.5)
                end
            end
        end
        
        X0tmp = X0mat;
        Tn = t(1:plot_int);
        
        % First run to get started
        for k = 1:nc
            subplot(h_s{k})
            hold on
            if plotstrt == 1
                X0 = X0mat(2:end,k);
                plot(X0(stx),X0(sty),[plotcol(k),'s'],'MarkerFaceColor',plotcol(k))
            end
            
            X = call_deri(Tn,par,X0mat(:,k),glo); % use call_deri.m to provide the output
            subplot(h_s{k})
            hold on
            hl1(k)=plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trainsz);
            if trailsz > 0
                plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trailsz);
            end
            drawnow
            X0tmp(2:end,k) = (X(end,:))';
            h_leg = legend(hl1(k),[leglab1,num2str(X0mat(1,k)),leglab2]); % create legend entries)
            set(h_leg,ft.name,ft.legend); % use standard font formats
        end
        
        % Continue from the first run for the remainder of the time vector
        for i = 2:floor(length(t)/plot_int) % note that some points might not be plotted now!
            for k = 1:nc
                Tn = t((i-1)*plot_int:i*plot_int); % new time vector
                X = call_deri(Tn,par,X0tmp(:,k),glo); % use call_deri.m to provide the output
                subplot(h_s{k})
                hold on
                delete(hl1(k)) % delete old plot segment before plotting a new one
                hl1(k)=plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trainsz);
                if trailsz > 0
                    plot(X(:,stx),X(:,sty),[plotcol(k),'-'],'LineWidth',trailsz);
                end
                drawnow
                X0tmp(2:end,k) = (X(end,:))'; % use previous last as new first
            end
        end
        if trailsz > 0
            for k = 1:nc
                delete(hl1(k)) % delete old plot segment at the end
            end
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
        
    case 7
        
        % =================================================================
        % make one 3d plot for each scenario
        % =================================================================
        for k = 1:nc % run through scenarios
            h_s{k} = subplot(m,n,k);
            hold on
            view(-20,25) % camera angle
            grid on
            xlabel(ylab(stx),ft.name,ft.label)
            ylabel(ylab{sty},ft.name,ft.label)
            zlabel(ylab{stz},ft.name,ft.label)
            % title(['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter', 'none',ft.name,ft.title)
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            
            if ~isempty(Xeq)
                for ieq = 1:size(Xeq,1)
                    stem3(Xeq(ieq,stx),Xeq(ieq,sty),Xeq(ieq,stz),'ko','MarkerFaceColor','y','BaseValue',axrange(stz,1),'LineWidth',2)
                end
            end
            
            if ~isempty(axrange) % if there is an axis range, use it
                axis([axrange(stx,:) axrange(sty,:) axrange(stz,:)])
            end
        end
        
        X0tmp = X0mat;
        Tn = t(1:plot_int);
        
        % First run to get started
        for k = 1:nc
            subplot(h_s{k})
            hold on
            if plotstrt == 1
                X0 = X0mat(2:end,k);
                plot3(X0(stx),X0(sty),X0(stz),[plotcol(k),'s'],'MarkerFaceColor',plotcol(k)) % plot starting point
            end
            
            X = call_deri(Tn,par,X0mat(:,k),glo); % use call_deri.m to provide the output
            subplot(h_s{k})
            hold on
            hl1(k)=plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trainsz);
            if trailsz > 0
                plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trailsz);
            end
            drawnow
            X0tmp(2:end,k) = (X(end,:))';
            h_leg = legend(hl1(k),[leglab1,num2str(X0mat(1,k)),leglab2]); % create legend entries)
            set(h_leg,ft.name,ft.legend); % use standard font formats
            set(h_leg,'Location','north')
        end
        
        % Continue from the first run for the remainder of the time vector
        for i = 2:floor(length(t)/plot_int) % note that some points might not be plotted now!
            for k = 1:nc
                Tn = t((i-1)*plot_int:i*plot_int); % new time vector
                X = call_deri(Tn,par,X0tmp(:,k),glo); % use call_deri.m to provide the output
                subplot(h_s{k})
                hold on
                delete(hl1(k)) % delete old plot segment before plotting a new one
                hl1(k) = plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trainsz);
                if trailsz > 0
                    plot3(X(:,stx),X(:,sty),X(:,stz),[plotcol(k),'-'],'LineWidth',trailsz);
                end
                drawnow
                X0tmp(2:end,k) = (X(end,:))'; % use previous last as new first
            end
        end
        if trailsz > 0
            for k = 1:nc
                delete(hl1(k)) % delete old plot segment at the end
            end
        end
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
end
