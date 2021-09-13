function [par_tmp,Xing,loglik_best] = auto_profiles(par,opt_prof,opt_optim)

% Usage: [par_tmp,Xing,loglik_best] = auto_profiles(par,opt_prof,opt_optim)
% 
% A helper function to run profiles for all fitted parameters, and
% re-optimise when locating a better value (returned in <par_tmp>).
%
% Author     : Tjalling Jager
% Date       : Dec. 2020
% Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 h_txt

% Set the number of sub-optimisations for a comparative analysis (maximum
% of three). This is only used when opt_prof.subopt = -1.
subopt_set = [0 10]; % compare zero and 10 sub-optimisations
% subopt_set = [0 20 40]; % different comparison (at this moment, 4 max!)

names   = glo2.names; % extract all field names of par
nfields = length(names);

drawnow % plot things that still need to be plotted
par_tmp = par; % work on a copy of the parameter structure
pmat    = packunpack(1,par_tmp,0);  % transform structure into a regular matrix
ind_fit = find(pmat(:,2)==1);
n_fit   = length(ind_fit);
loglik_best = -1; % this will collect better minloglik, if found

% Reset options for the profiling specifically for this automated analysis
opt_prof.brkprof  = 1; % set to 1 to break the profiling when a better optimum is located
opt_prof.verbose  = 0; % set to 0 to suppress all output from calc_proflik to screen

switch opt_prof.subopt
    case -1 % do the comparison between different number of sub-optimisations
        comp_sub = 1;
    case -2 % do comparison with simulated annealing sub-optimisations
        comp_sub = 2;
        subopt_set = [0 1];
    otherwise
        comp_sub = 0;
        subopt_set = opt_prof.subopt;
end

disp(' ')
disp('Starting automatic calculations.')
if opt_prof.verbose == 0
    disp('No results will be printed on screen or plotted until the analysis is finished.')
end

%% Loop profiling and optimisation until no better value results

prof_coll = cell(length(subopt_set),nfields);

for jj = 1:length(subopt_set) % run through the set of sub-optimisations
    
    switch comp_sub
        case 1
            opt_prof.subopt = subopt_set(jj); % take jj-th element of sub-optimisation number
            disp(['Starting next run with ',num2str(subopt_set(jj)),' sub-optimisations'])
        case 2
            opt_prof.subann = subopt_set(jj); % set to 1 to use simulated annealing (followed by simplex) instead of suboptimisations
            if jj == 1
                disp('Starting next run without sub-optimisations')
            else
                disp('Starting next run with sub-optimisations using annealing')
            end
    end
    
    order  = ind_fit; % start with the profiles in order of the par structure
    hurrah = 0; % flag that we can stop the analysis and print the results
    
    while hurrah == 0
        
        % Prepare/empty cell arrays to catch the profile output
        Xing = cell(1,nfields);
        
        for ii = 1:n_fit % run through all parameters
            i = order(ii); % this helps to modify the order later if needed
            disp(['  Starting a profile for parameter ',names{i},' (',num2str(ii),' of ',num2str(n_fit),')'])
            [a,par_better,d,~] = calc_proflik(par_tmp,names{i},opt_prof); % calculate a profile
            if isstruct(par_better) % profiles have not finished and returned a better par structure
                hurrah = 0; % we cannot stop just yet ... better optimum is located
                break % break from the for loop
            else
                Xing{i} = a; % collect the intervals coming out
                prof_coll{jj,i} = d; % collect the profiles coming out
            end
            hurrah = 1; % no better optimum, so move to next parameter or end
        end
        
        if hurrah == 0 % the profiling has stopped as a better value was found
            opt_optim.it = 0; % show iterations of the simplex optimisation (1, default) or not (0)
            warning('off','backtrace')
            warning('Better optimum was found, initiating new optimisation')
            disp(' '), warning('on','backtrace')
            [par_tmp,loglik_best] = calc_optim(par_better,opt_optim); % new par_tmp
            order(order==i) = []; % remove the parameter from the list ...
            order = [i;order]; % and put it in first place so it is the first one profiled again
        end
        
    end
       
end

%% Display the results on screen and to the diary
% Note: if comparing with and without sub-opts, the run with
% sub-optimisations is printed!

% ready to display results
diary (glo.diary) % collect screen output in the diary "results.out"
disp(' ')
if comp_sub > 0
    disp(['Parameters and bounds are for the last runs with ',num2str(subopt_set(end)),' sub-optimisations.'])
end
disp('Parameter, best fit value, interval')
disp('==========================================================')
pmat = packunpack(1,par_tmp,0);  % transform structure into a regular matrix
for i = 1:nfields % run through all parameters
    if isempty(Xing{i}) % then it was not a fitted parameter
        fprintf('%-6s %10.4g (parameter not fitted) \n',names{i},pmat(i,1))
    else
        % column 3 is for warnings when hitting a boundary
        if Xing{i}(1,2) == 1 % interval open on lower end
            disp('Warning: taking the lowest parameter value as lower confidence limit')
        end
        if Xing{i}(end,2) == 1 % interval open on upper end
            disp('Warning: taking the highest parameter value as upper confidence limit')
        end
        if size(Xing{i},1) == 2 % no problems, it is a single interval
            fprintf('%-10s best: %#10.4g interval: %#10.4g - %#1.4g \n',names{i},pmat(i,1),Xing{i}(1,1),Xing{i}(2,1))
        else
            disp('The confidence interval is a broken set (check likelihood profile to check these figures)')
            fprintf('%-10s best: %#10.4g interval: \n',names{i},pmat(i,1))
            for ix = 1:size(Xing{i},1)/2
                fprintf('   interval %1.1g: %#10.4g - %#1.4g \n',ix,Xing{i}((ix-1)*2+1,1),Xing{i}((ix-1)*2+2,1))
            end
        end
    end
end
disp('==========================================================')
disp(['Time required: ' secs2hms(toc)])
disp(' ')
diary off  % close results.out

%% Make a multiplot with all final profiles

% Calculate size of multiplot
n = ceil(sqrt(n_fit));
m = ceil(n_fit/n);
[figh,ft] = make_fig(m,n); % create a figure window of correct size

chicrit = 3.8415; % = chi2inv(0.95,1); % the chi square criterion for 1 df and 95%

% for plotting, put pmat on log-scale where needed
pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));

switch length(subopt_set)
    case 1
        line_set = {'k.-'};
    case 2
        line_set = {'k.:','k.-'};
    case 3
        line_set = {'k.:','k.--','k.-'};
    case 4
        line_set = {'k.:','k.--','k.-.','k.-'};
end

g_rem = cell(n_fit,1);
for i = 1:n_fit
    g = subplot(m,n,i); % subplot 
    hold on
    h1 = gca; % remember the current axis number for plotting!
    set(h1,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
    g_rem{i} = g; % remember subplot handle
    
    if pmat(ind_fit(i),5) == 0
        xlabel(h1,['parameter ',names{ind_fit(i)},' (log-scale)'],ft.name,ft.label)
    else
        xlabel(h1,['parameter ',names{ind_fit(i)}],ft.name,ft.label)
    end
    ylabel(h1,'minus 2x log-likelihood ratio',ft.name,ft.label)
    
    % And make a nice plot
    for jj = 1:length(subopt_set) % this is for the comparison between 0 and 10 sub-opts
        h_line{jj} = plot(h1,prof_coll{jj,ind_fit(i)}(:,1),prof_coll{jj,ind_fit(i)}(:,2),line_set{jj});
        Le{jj} = [num2str(subopt_set(jj)),' sub-opt.'];
    end
    
    % Plot the chi-square criterion for 95% and 1 degree of freedom
    plot(h1,[min(prof_coll{end,ind_fit(i)}(:,1)) max(prof_coll{end,ind_fit(i)}(:,1))],[chicrit chicrit],'k:')
    plot(h1,pmat(ind_fit(i)),0,'ko','MarkerFaceColor','w') % plot the max lik value as a circle
    ylim([0 11]) % limit y-axis to just above the cut-off criterion
    
    a = [min(prof_coll{end,ind_fit(i)}(:,1)) max(prof_coll{end,ind_fit(i)}(:,1))];
    a = [a(1)-diff(a)/20 a(2)+diff(a)/20];
    xlim(a) % limit x-axis, based on last run
end

if comp_sub > 0 % then we also have an array with legend entries
    h_leg = legend(Le);
    set(h_leg,ft.name,ft.legend);
end

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
h_txt = text(0.5, 1,['Plotted from: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
set(h_txt,ft.name,ft.text); % use standard formatting for this header

if comp_sub > 0 && length(subopt_set) == 2
    for i = 1:n_fit % compare the two profiles
        if sum(size(prof_coll{2,ind_fit(i)}) == size(prof_coll{1,ind_fit(i)}))==2 % they are of equal size
            diff_prof_x = prof_coll{2,ind_fit(i)}(:,1)./prof_coll{1,ind_fit(i)}(:,1); % relative difference
            diff_prof_y = prof_coll{1,ind_fit(i)}(:,2)-prof_coll{2,ind_fit(i)}(:,2); % absolute difference
            if max(diff_prof_x) < 1.01 && min(diff_prof_x) > 0.99 % the x-axis is close enough
                if max(abs(diff_prof_y)) < 0.1 % the y-axis is probably close enough
                    fprintf('very similar profile for parameter %-4s (max. diff. %1.2g) \n',names{ind_fit(i)},max(abs(diff_prof_y)))
                else
                    fprintf('profiles seem to differ for parameter %-4s (y-values too different, max. diff. %1.2g) \n',names{ind_fit(i)},max(abs(diff_prof_y)))
                end
            else
                fprintf('profiles seem to differ for parameter %-4s (x-values too different) \n',names{ind_fit(i)})
            end
        else
            fprintf('profiles seem to differ for parameter %-4s (no. points differs) \n',names{ind_fit(i)})
        end
    end
    disp(' ')
end

%% Save plot if needed
    
if glo.saveplt > 0 % if we want to save the plot
    savenm = ['prof_all_',glo.basenm];%
    save_plot(figh,savenm,h_txt);
end

