function calc_localsens(par_plot,opt_sens)

% Usage: calc_localsens(par_plot,opt_sens)
% 
% This function performs a local sensitivity analysis of the model. All
% model parameters are increased one-by-one by a small percentage. The
% sensitivity score is either the scaled version (dX/X p/dp) or the
% absolute version (dX p/dp).
%
% Options to set in a structure <opt_sens>, which is filled in
% <prelim_checks> and can be overwritten in main script, where needed.
%
% Author     : Tjalling Jager 
% Date       : May 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 X0mat

% extract parameters from the general globals, glo and glo2
t     = glo.t;
n_X   = glo2.n_X;
names = glo2.names;
n_s   = size(X0mat,2); % number of scenarios

% read options from structure
sens_type = opt_sens.type;
Fchange   = opt_sens.step;
if opt_sens.state(1) == 0
    state_sel = 1:n_X; % use all state variables for the sensitivity analysis
else
    state_sel = opt_sens.state;
end
savenm1 = glo.basenm; % scriptname for title on sensitivity analyses

pmat = packunpack(1,par_plot,0); % use the fitted parameters

% % Now, I use ALL parameters in the analyses, also the ones that were not fitted
% % For future modification, I keep these lines below in comments.
% pfit    = pmat(pmat(:,2)==1,1); % parameter values that are to be estimated
% ind_fit = find(pmat(:,2)==1);   % index to fitted parameters

%% Start by calculating the standard model output (with best-fit parameters)

% initialise the matrices that will collect the basic output (times x scenarios)
Xbase = cell(n_X,1); % pre-define structure
for i = 1:n_X % run through state variables
    Xbase{i} = zeros(length(t),n_s); % initialise state variable output with zeros
end
for j = 1:n_s % run through our scenarios
    Xout = call_deri(t,par_plot,X0mat(:,j)); % use call_deri.m to provide the output for one scenario
    for i = 1:n_X  % run through state variables
        Xbase{i}(:,j) = Xout(:,i); % collect state variable i into structure Xbase
    end
end

%% Calculate the model output with slightly changed parameters, and the sensitivty coefficients

% initialise the matrices that will collect the output (times x scenarios x parameter counter)
sens = cell(n_X,1); % pre-define structure
for i = 1:n_X % run through state variables
    sens{i} = zeros(length(t),n_s,size(pmat,1)); % initialise sensitivities matrix with zeros
end
for p = 1:size(pmat,1) % run through all parameters in the script
    pmat_tmp = pmat; % start from fresh parameter matrix
    pmat_tmp(p,1) = pmat_tmp(p,1) * (1+Fchange); % change one parameter slightly
    par_k = packunpack(2,0,pmat_tmp); % transform parameter matrix into a structure
    for j = 1:n_s % run through our scenarios
        Xout = call_deri(t,par_k,X0mat(:,j)); % use call_deri.m to provide the output for one scenario
        for i = 1:n_X  % run through all state variables
            switch sens_type
                case 1
                    relstate = (Xout(:,i) - Xbase{i}(:,j)) ./ Xbase{i}(:,j); % relative change in state
                case 2
                    relstate = (Xout(:,i) - Xbase{i}(:,j)); % absolute change in state
            end
            sens{i}(:,j,p) = relstate / Fchange; % sensitivity coefficient
        end
    end
end
      
%% And finally, plot the results in a multiplot (one panel for each scenario)

for ii = 1:length(state_sel) % make a new figure window for each state variable
    % only do this for the states we selected at the start of this script
    
    plotcol = 'krbgcm'; % 6 colors should be enough ...
    % make one plot with subplots for each scenario (i.e., concentration)
    n = ceil(sqrt(n_s));
    m = ceil(n_s/n);
    
    [figh,ft] = make_fig(m,n); % create figure of correct size
    
    sens_max = max(max(max(sens{state_sel(ii)})));
    sens_min = min(min(min(sens{state_sel(ii)})));
    
    for j = 1:n_s % run through scenarios
        if n_s > 2 % for 1 or 2 scenarios, making the plot larger is no fun
            g = subplot(m,n,j);
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
            p = get(g,'position');
            p([3 4]) = p([3 4])*1.10; % add 10 percent to width and height
            set(g, 'position', p);
        else
            set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        end
        if n_s > 1
            title([glo.leglab1,num2str(X0mat(1,j)),' ',glo.leglab2],ft.name,ft.title)
        end
            
        hold on
        
        for p = 1:size(pmat,1) % run through all parameters in the script
            if p <= 6
                plot(t,sens{state_sel(ii)}(:,j,p),[plotcol(p),'-'],'LineWidth',2)
            else % for more than 6 parameters, use dotted line
                plot(t,sens{state_sel(ii)}(:,j,p),[plotcol(p-6),':'],'LineWidth',2)
            end
        end
        ax = gca; % handle to current axis system
        if j>(n*(m-1)) % only put xlabel on bottom row
            xlabel(glo.xlab,ft.name,ft.label)
        else
            % ax.XTickLabel = {[]}; % this is okay for new Matlab versions
            set(ax,'XTickLabel',[]); % this works with old versions as well
        end
        if j==1 || (j-1)/n == floor((j-1)/n) % only put y labels on first column
            switch sens_type
                case 1
                    ylabel(['relative sens. for state ',num2str(state_sel(ii))],ft.name,ft.label)
                case 2
                    ylabel(['absolute sens. for state ',num2str(state_sel(ii))],ft.name,ft.label)
            end
        else
            % ax.YTickLabel = {[]};
            set(ax,'YTickLabel',[]);
        end
        
        plot([t(1) t(end)],[0 0 ],'k:','LineWidth',2)
        axis([t(1) t(end) sens_min-0.01*(sens_max-sens_min) sens_max]) % set all axis within the multiplot the same
    end
    legend(names)
    
    if opt_sens.notitle == 0
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Local sensitivities for state ',num2str(state_sel(ii))],'HorizontalAlignment','center','VerticalAlignment', 'top');
        set(h_txt,ft.name,ft.text); % use standard formatting for this header
    end
    
    if glo.saveplt > 0 % if we want to save the plot
        savenm = ['localsens_',savenm1,'_state_',num2str(ii)];%
        save_plot(figh,savenm,h_txt);
    end
end
