function figh = plot_grid(pmat,coll_ok,coll_all,coll_prof_pruned,figh,SETTINGS_OPTIM)

% Usage: figh = plot_grid(pmat,coll_ok,coll_all,coll_prof_pruned,figh,SETTINGS_OPTIM)
%
% Function to make a multipanel plot of the sample as generated in
% <calc_parspace>. The panels on the diagonal are for the profile
% likelihoods (for single parameters). The panels under the diagonal are
% for each combination of two parameters. This function is called from
% <calc_parspace> and <calc_optim_ps>. From <calc_parspace>, it will be called
% several times as the optimisation progresses (unless the option
% <opt_optim.ps_plots> is set to zero).
%
% Inputs
% <pmat>        parameter matrix
% <coll_ok>     part of the sample that will be taken to the next round (matrix)
% <coll_all>    the entire sample (matrix)
% <coll_prof_pruned> information (cell array) to make the profile likelihoods (red line on diagonal plots)
% <figh>        handle to the figure, to ensure that updates can be plotted in the same figure
% <SETTINGS_OPTIM> = settings read from setup_settings by <calc_optim_ps>
%
% Outputs
% <figh>        handle to the figure, to ensure that updates can be plotted in the same figure
% 
% Author     : Tjalling Jager 
% Date       : September 2020
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <plot_grid> code that is
% distributed as part of the Matlab version of openGUTS (see
% http://www.openguts.info). Therefore, this code is distributed under
% the same license as openGUTS (GPLv3). The modifications are not in the
% algorithm itself but only to ensure that the code operates in the general
% BYOM framework.
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%  
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

%% BLOCK 1. Initial things

global glo2

names = glo2.names; % names of the model parameters

% BLOCK 1.1. Extract useful information from <pmat>.
ind_fit  = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
n_fit    = length(ind_fit);      % number of fitted parameters
ind_log  = find(pmat(:,5) == 0); % indices to log-scale parameters in total <pmat> (vector)
bnds_tmp = pmat(:,[3 4]);        % copy parameter bounds (matrix)
bnds_tmp(ind_log,[1 2]) = log10(bnds_tmp(ind_log,[1 2])); % and also log-transform bounds where needed
bnds_tmp = bnds_tmp(ind_fit,:);  % only keep the ones for the fitted parameters

% BLOCK 1.2. Chi2-criteria and indices to parameter sets to plot for inner and outer rim.
% We work here with the log-likelihood itself, so the chi2 criterion needs to be divided by 2.
chicrit_joint  = 0.5 * SETTINGS_OPTIM.crit_table(n_fit,1); % criterion for joint 95% CI or parameters
chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1);     % criterion for single-parameter CIs

% Find indices to sets within inner and outer rim.
coll_best  = coll_all(1,:); % this is the best-fitting parameter set (<coll_all> is sorted)
ind_single = find(coll_all(:,end) < coll_best(end) + chicrit_single,1,'last'); % index to last element of sample that is still in inner rim
ind_fin95  = find(coll_all(:,end) < coll_best(end) + chicrit_joint,1,'last');  % index to last element of sample that is still in joint CI

% The <coll_ok> is only plotted during the initial rounds of the
% optimisation, to provide a sense of the progress. For the last round, it
% is not input to this function, and we will plot the joint 95% set
% instead.
if isempty(coll_ok) % if there is no <coll_ok> provided ...
    coll_ok = coll_all(ind_single+1:ind_fin95,:); % use the outer rim for those points (the ghost symbols)
    % Note: I now exclude the inner rim, to avoid plotting them twice.
end
coll_inner = coll_all(1:ind_single,:); % inner rim for plotting

%% BLOCK 2. Prepare figure window and calculate plotting bounds
% Some care is needed in this multi-panel plot to make sure that all panels
% witin a column have the same x-axis scaling, and that all binary plots on
% the same row have the same y-axis.

if ~isempty(figh) % if a handle is provided in the input ...
    figure(figh)  % use existing handle and make sure that the multiplot is the current plot
    clf(figh)     % clear anything that was plotted previously
    [~,ft] = make_fig(1,2,1); % call make_fig just to return the font properties
else % otherwise ...
    if n_fit > 1
        [figh,ft] = make_fig(n_fit,n_fit); % make empty figure window of correct size
    else % if only one parameter is fitted ...
        [figh,ft] = make_fig(1,2); % make figure window of 1x2 (this is treated differently)
    end
end
hold on % keep plotted info even when there are subsequent plotting events

% First caculate the axis bounds we will use for the axes with parameters.
% We use <coll_ok> for that, which is the largest set we'll print (either
% the sets to continue to the next round or the joint CI). Note that last
% column in <coll_ok> is the minloglik, which we can exclude. Also note
% that <min_sample> and <max_sample> will be vectors (one value for each
% parameter): the 1 in the calls signals that <min> and <max> need to
% operate over all rows per column (which is Matlab default, so
% superfluous).
min_sample1 = min(coll_ok(:,1:end-1),[],1); % (row vector)
max_sample1 = max(coll_ok(:,1:end-1),[],1); % (row vector)

% However, in some extreme cases, the best value (and/or part of inner rim)
% may be outside <coll_ok> (when the best value has jumped to another
% location). Therefore, it is good to also see of the inner rim is not
% extending beyond the min/max bounds calculated.
min_sample2 = min(coll_inner(:,1:end-1),[],1); % (row vector)
max_sample2 = max(coll_inner(:,1:end-1),[],1); % (row vector)

% Take the outer edges of the two ranges calculated.
min_sample = min([min_sample1;min_sample2],[],1); % (row vector)
max_sample = max([max_sample1;max_sample2],[],1); % (row vector)

% We need to adjust these bounds for plotting profiles. Profile likelihoods
% may extend beyond the range of the sample, and we need to see that.
if ~isempty(coll_prof_pruned) % when there is profile data ...
    for i_p = 1:n_fit % run through the fitted parameters
        % Adjust the bounds based on the info in the profile.
        if ~isempty(coll_prof_pruned{i_p}) % NEW! In extreme cases it can be empty!
            min_sample(i_p) = min(min_sample(i_p),min(coll_prof_pruned{i_p}(:,1)));
            max_sample(i_p) = max(max_sample(i_p),max(coll_prof_pruned{i_p}(:,1)));
        end
    end
end

% Add a little bit extra for the plotting bounds so that the plots look
% nice and not cluttered on top of the axes.
diff_sample = max(abs(0.01*min_sample),abs(max_sample-min_sample)); % calculate the total range of the sample
% I cannot remember why I want a minimum range of abs(0.01*min_sample) ...
% this is probably only needed for extreme cases.
min_sample  = min_sample - 0.05 * diff_sample; % lower the minimum range by x%
max_sample  = max_sample + 0.05 * diff_sample; % increase the maximum range by x%

%% BLOCK 3. Plot two-dimensional projections of the parameter sample.

if n_fit > 1 % only when fitting more than two parameters (for one fitted parameter, a special plot is made)
    for i_1 = 1:n_fit-1 % this will run through columns in the figure
        for i_2 = i_1+1:n_fit % this will run through rows in figure (for the panels below the diagonal)
            
            g  = subplot(n_fit,n_fit,(i_2-1)*n_fit+i_1); % construct a sub-plot at the correct panel
            ax = gca; % handle to the current axis system
            set(ax,'LineWidth',1,ft.name,ft.ticks) % set axes to the standard font format
            
            % Modify plot to decrease white space between panels (looks better in Matlab).
            p = get(g,'position'); % get the size of the plot
            p([3 4]) = p([3 4])*1.25; % add some fraction to width and height
            if n_fit == 2 % when plotting just 2 pars, need to move them down a bit ...
                p(2) = p(2)-0.04;
            end
            set(g, 'position', p); % set the plot to the new size
            hold on % make sure that hold is on, so we can plot multiple things
            
            % Plot the sample points in there: inner rim marked in green.
            plot(coll_ok(:,i_1),coll_ok(:,i_2),'ko','MarkerFaceColor','w','MarkerEdgeColor','c') % plot all tries to continue with (or outer rim) as cyan open circles
            plot(coll_inner(:,i_1),coll_inner(:,i_2),'ko','MarkerFaceColor','g') % plot the ones from inner rim in green
            
            % Plot bounds on the parameters as blue lines.
            plot([bnds_tmp(i_1,1) bnds_tmp(i_1,1)],[bnds_tmp(i_2,1) bnds_tmp(i_2,2)],'b-','LineWidth',1.5)
            plot([bnds_tmp(i_1,2) bnds_tmp(i_1,2)],[bnds_tmp(i_2,1) bnds_tmp(i_2,2)],'b-','LineWidth',1.5)
            plot([bnds_tmp(i_1,1) bnds_tmp(i_1,2)],[bnds_tmp(i_2,1) bnds_tmp(i_2,1)],'b-','LineWidth',1.5)
            plot([bnds_tmp(i_1,1) bnds_tmp(i_1,2)],[bnds_tmp(i_2,2) bnds_tmp(i_2,2)],'b-','LineWidth',1.5)
            
            % Plot best parameter set in yellow.
            plot(coll_best(i_1),coll_best(i_2),'ko','MarkerFaceColor','y','MarkerSize',8) % plot the best fit parameter set so far in yellow
            
            % Set axis bounds to the overal min-max for these parameters.
            axis([min_sample(i_1) max_sample(i_1) min_sample(i_2) max_sample(i_2)]) % set the axis bounds
            
            % Figure out when to put labels on the axis (only outer
            % panels). Labels are made with standard font formatting.
            if i_1 == 1 % we are in the first column
                if pmat(ind_fit(i_2),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                    ylabel([names{ind_fit(i_2)},' (log)'],ft.name,ft.label) % so indicate log in the axis label
                else
                    ylabel([names{ind_fit(i_2)}],ft.name,ft.label)
                end
                if i_2 < n_fit % if we are not on the last row, remove axis numbers
                    set(ax,'XTickLabel',[]); 
                else % otherwise add an x-axis label
                    if pmat(ind_fit(i_1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                        xlabel([names{ind_fit(i_1)},' (log)'],ft.name,ft.label) % so indicate log in the axis label
                    else
                        xlabel([names{ind_fit(i_1)}],ft.name,ft.label)
                    end
                end
            elseif i_2 == n_fit % we are in the last row
                if pmat(ind_fit(i_1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                    xlabel([names{ind_fit(i_1)},' (log)'],ft.name,ft.label) % so indicate log in the axis label
                else
                    xlabel([names{ind_fit(i_1)}],ft.name,ft.label)
                end
                if i_1>1 % if we are not on the first column, remove axis numbers
                    set(ax,'YTickLabel',[]); 
                end
            else % otherwise remove numbers on both axes
                set(ax,'YTickLabel',[]); 
                set(ax,'XTickLabel',[]); 
            end
            
        end
    end
end

%% BLOCK 4. Create a legend in top-right sub-plot.
% This requires a bit of trickery to get right in Matlab.

if n_fit > 1 % if we have more than 1 fitted parameter ...
    h_annot = subplot(n_fit,n_fit,n_fit); % make an extra subplot, and remember the handle
    % Modify plot to decrease white space between the panels.
    p = get(h_annot,'position');
    p([3 4]) = p([3 4])*1.2; % add to width and height
    if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
        p(2) = p(2)-0.04;
    end
    set(h_annot, 'position', p); % this minimises the white space between the panels
else % if we have one fitted parameters, the legend will be next to the profile plot
    h_annot = subplot(1,2,2); % make an extra subplot, and remember the handle
end
h1  = gca; % remember the current axis
hold on

% Plot fake data points, with same symbols as in regular plots (to put them
% into the legend automatically).
if n_fit > 1
    plot(0,0,'ko','MarkerFaceColor','w','MarkerEdgeColor','c') % plot all tries to continue with as ghost circles
    plot(0,0,'ko','MarkerFaceColor','g') % plot the ones from inner rim in green
end
plot(0,0,'ko','MarkerFaceColor','y','MarkerSize',8) % plot the best fit parameter set so far in yellow
plot([0 0],[0.1 0.1],'b-','LineWidth',1.5)

% Define text for legend entries.
if isempty(coll_prof_pruned) && n_fit > 1 % then we are in an intermediate round
    legend_txt = {'to next round','inner rim','best fit set','space bounds'};
else % them we are in a final round, so also plot symbols for the pseudo-profiles
    plot(0,0,'c.') % all values within joint 95% CI (df = fitted pars)
    plot(0,0,'g.') % all values within single par CI (df = 1)
    plot([0 0],[0.1 0.1],'k-') % line for edges of single par CI
    plot([0 0],[0.1 0.1],'r-','LineWidth',1.5)
    if n_fit > 1
        legend_txt = {'joint 95% CI','inner rim (df=1)', ...
            'best fit set','space bounds','joint 95% CI', ...
            'inner rim (df=1)','cut-off (df=1)','profile refine'};
        
        if n_fit > 5 && chicrit_joint < 6 % if we limited the outer rim ...
            legend_txt{1} = 'outer rim (df=5)'; % also modify the legend entry
            legend_txt{5} = 'outer rim (df=5)'; % also modify the legend entry
        end
            
    else
        legend_txt = {'best fit set','space boundaries','joint 95% CI','inner rim (df=1)' ...
            'cut-off for CIs','profile refinement'};
    end     
end
h_leg_tot = legend(legend_txt); % create a legend with all entries
set(h_leg_tot,ft.name,ft.legend); % use standard font formats

xlim([1 10]) % make sure the 'data' are off screen, so invisible
dimleg = h_leg_tot.Position; % extract position of legend made
dim    = h_annot.Position; % the position of the sub-plot is also the correct place for the legend
% Move the legend to the new sub-plot and set it to top-left of panel.
h_leg_tot.Position = [dim(1) dim(2)+dim(4)-dimleg(4) dimleg(3:4)];
h1.Position = [dim(1) dim(2)+dim(4)-dimleg(4) 0.01 0.01];
h1.Visible  = 'off'; % hide axes

%% BLOCK 5. Plot the profile likelihoods on the diagonal.

if ~isempty(coll_prof_pruned) || n_fit == 1 % only when <coll_prof> is provided, or when there is only 1 fitted parameter
    % We also do that if there is only one fitted parameter. For that case,
    % we skipped the binary parameter plot as that makes no sense.
    
    for i_p = 1:n_fit % run through all fitted parameters
        
        figure(figh) % make sure that the multiplot is the current plot
        if n_fit > 1 % for more than 1 fitted parameter ...
            g = subplot(n_fit,n_fit,(i_p-1)*n_fit+i_p); % create an empty sub-plot on diagonal
            % Modify plot to decrease white space between panels.
            p = get(g,'position'); % get the size of the plot panel
            p([3 4]) = p([3 4])*1.25; % add some percentage to width and height
            if n_fit == 2 % when plotting just 2 pars, need to move them down a bit
                p(2) = p(2)-0.04;
            end
            set(g, 'position', p); % set plot to new size
        else % for one fitted parameter, we have a 2-panel plot, so take first sub-plot
            g = subplot(1,2,1); % sub-plot on first position
        end
        h1 = gca; % remember the current axis number for plotting!
        set(h1,'LineWidth',1,ft.name,ft.ticks) % adapt to standard axis formatting of fonts
        
        % Determine which plots get a label on which axis.
        parnum  = ind_fit(i_p); % find the index for fitted parameter <i_p>
        if i_p == 1 % if we have the first fitted parameter ...
            ylabel(h1,'loglik ratio',ft.name,ft.label) % put a label on the y-axis
            if n_fit > 1 % for more than 1 fitted parameter ...
                set(h1,'XTickLabel',[]); % remove tick labels on x-axis
            else % for 1 parameter, label the x-axis with the parameter name
                if pmat(parnum,5) == 0
                    xlabel(h1,[names{parnum},' (log)'],ft.name,ft.label)
                else
                    xlabel(h1,[names{parnum}],ft.name,ft.label)
                end
            end
        elseif i_p == n_fit % if we have reached the last plot, we need an x-axis label
            if pmat(parnum,5) == 0
                xlabel(h1,[names{parnum},' (log)'],ft.name,ft.label)
            else
                xlabel(h1,[names{parnum}],ft.name,ft.label)
            end
            set(h1,'YTickLabel',[]); % remove tick labels on y-axis
        else
            set(h1,'YTickLabel',[]); % remove tick labels on y-axis
            set(h1,'XTickLabel',[]); % remove tick labels on x-axis
        end
        
        hold on
        
        % Plot the sample itself as points (log-likelihood ratio versus parameter value); same colours as in binary plots.
        plot(coll_ok(:,i_p),coll_ok(:,end)-coll_best(end),'c.') % all values within joint CI (excl. inner rim)
        plot(coll_inner(:,i_p),coll_inner(:,end)-coll_best(end),'g.') % all values within inner rim

        % Plot some horizontal lines to indicate the edge of the CI and the
        % band that will be propagated to model predictions.
        plot([min_sample(i_p) max_sample(i_p)],[chicrit_single chicrit_single],'k-') % line for edges of single par CI
        plot([min_sample(i_p) max_sample(i_p)],0.5*[SETTINGS_OPTIM.crit_prop(1) SETTINGS_OPTIM.crit_prop(1)],'k:') % lower edge of propagation interval
        plot([min_sample(i_p) max_sample(i_p)],0.5*[SETTINGS_OPTIM.crit_prop(2) SETTINGS_OPTIM.crit_prop(2)],'k:') % upper edge of propagation interval
         
        % Plot the profile refinement as red line.
        if ~isempty(coll_prof_pruned) % this skips plotting a profile refinement for the 1-fitted-parameter case
            plot(coll_prof_pruned{i_p}(:,1),coll_prof_pruned{i_p}(:,2)-coll_best(end),'r-','LineWidth',1.5) % plot the profile refinment as red line
        end
        
        % Plot bounds on the parameter space as blue lines.
        plot([bnds_tmp(i_p,1) bnds_tmp(i_p,1)],[0 chicrit_joint],'b-','LineWidth',1.5)
        plot([bnds_tmp(i_p,2) bnds_tmp(i_p,2)],[0 chicrit_joint],'b-','LineWidth',1.5)
        
        % Plot best parameter value with a yellow marker.
        plot(coll_best(i_p),0,'ko','MarkerFaceColor','y','MarkerSize',8) % best fit parameter value as yellow circle
        
        % Set the x-axis to the overal bounds, and the y-axis to the chi2 criterion.
        xlim([min_sample(i_p) max_sample(i_p)]);
        if n_fit > 1
            ylim([0 chicrit_joint]);
        else
            ylim([0 chicrit_joint+0.5]); % for 1 fitted parameter, take a bit extra
        end
            
        drawnow % update the plot
        
    end
    
end

% pause % handy if we want to look at each step in detail

drawnow % update the plot
