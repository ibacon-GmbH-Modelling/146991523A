function [par_best,XingS] = calc_likregion(par,nr_lhs,opt_likreg,opt_prof,varargin)

% Usage: [par_best,XingS] = calc_likregion(par,nr_lhs,opt_likreg,opt_prof,varargin)
%
% Calculates a likelihood-based joint confidence region. First, profile
% likelihoods for all fitted parameters are made to find the edges of the
% hyperbox that contains the true confidence region (of which the shape is
% then still unknown). Next, Latin Hypercube sampling is used to sample the
% hyperbox (shooting), and a likelihood-ratio test is used to decide which
% ones belong to the confidence region, and which are outside.
%
% This function also calculates the confidence intervals for the single
% parameters; it uses <calc_proflik.m> for that. It has to do the profiling
% anyway to get the edges of the 'hyperbox', so you get the intervals on
% the single parameters for free. Note that if you used <calc_proflik.m>
% already for all fitted parameters, you can used the saved sample
% (filenames with ..._LR.mat) and skip profiling for this function.
%
% Second input is the number of samples that we aim for in the inner rim
% (within the chi2 criterion for df=1). This will be the set that is to be
% used for creating CIs on model predictions. Note that the we add a little
% bit to the chi2 criterion: since we use a discrete sample to approximate
% a 'hyper hull' in parameter space.
% 
% This function should run without the statistics toolbox of Matlab. But
% the LHS sampling will then be replaced by normal random sampling.
%
% For possible options to set in a structure as third argument
% (<opt_likreg>) see <prelim_checks.m>. This function also needs to the
% options structure <opt_prof> for the profiling part, and optionally
% <opt_optim> to allow for refitting when a better optimum is found during
% profiling.
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 h_txt

% read options from structure
skip_prof = opt_likreg.skipprof; % skip profiling step; use boundaries from saved likreg set (1) or profiling (2)
chull     = opt_likreg.chull; % set to 1 to plot convex hull that approximates 95% edges of the likelihood region
axbnds    = opt_likreg.axbnds; % bind axes on the bounds of the hyperbox (1), accepted sample (2), or inner region (3)
burst     = opt_likreg.burst; % number of random samples from parameter space taken every iteration
lim_out   = opt_likreg.lim_out; % set to 1 to sample from a smaller part of space (enough for forward predictions)
brkprof   = opt_prof.brkprof; % set to 1 to stop profiling when better optimum is located, 2 to re-fit

if ~isempty(varargin) && ~isempty(varargin{1})
    opt_optim = varargin{1}; 
elseif brkprof == 2
    brkprof = 1; % we cannot refit without opt_optim
    if verbose ~= 0
        warning('off','backtrace')
        warning('No re-fitting will be done as opt_optim is left empty as input. Will break instead')
        warning('on','backtrace')
    end
end

%% Initial things

filenm      = glo.basenm;
par_best    = []; % make sure the output is defined, even if we stop prematurely
names       = glo2.names;

% By default, this file will use the name of the base script to load a MAT
% file. However, we also want to use saved MAT files for predictions in new
% script files. The glo.mat_nm can be used for that.
if isfield(glo,'mat_nm')
    filenm_load = glo.mat_nm;
else
    filenm_load = glo.basenm;
end

% list of all critical values of the chi-square for 95% with df 1 to 5
% that way we can work without the statistics toolbox in most cases
chitable = [3.8415 5.9915 7.8147 9.4877 11.07];

diary (glo.diary) % collect output in the diary "results.out"

disp(' ')
if exist('lhsdesign','file')~=2 % when lhsdesign does not exists as an m-file in the path
    warning('off','backtrace')
    Warning('You cannot use Latin-hypercube sampling; you need the statistics toolbox for that.')
    Warning('Instead, uniform random sampling will be used (more samples might be needed to obtain good coverage).')
    disp(' '), warning('on','backtrace')
end
if ~isfield(par,'tag_fitted') % apparently, parameters have not been fitted
    warning('off','backtrace')
    warning('Parameters have not been fitted, so results may not be very meaningful!')
    disp(' '), warning('on','backtrace')
end
    
%% Load info from saved MAT file if needed
% And check if it dffers from the par that is entered as input.

names_tmp = names; % work with a copy of names, as the saved set may have different names!
% if the saved set has different names, this likely will produce an error
% elsewhere anyway.
pmat_tmp  = packunpack(1,par,0); % use this to compare to the saved one, if re-using a saved set

if nr_lhs == -1 % do NOT make a new sample, but use the saved one!
    if exist([filenm_load,'_LR.mat'],'file') ~= 2
        error('There is no likelihood-region sample saved, so run calc_likregion again with a positive number of samples')
    end
    load([filenm_load,'_LR']) % load the random sample from the last likreg run
    % this loads the random sample in rnd, parameter matrix par and selection matrix par_sel
    % it also loads information from the profiling
    acc = rnd(rnd(:,end)<chitable(1),:); % put the sets in inner rim (df=1) into acc
    
    names_tmp  = fieldnames(par); % extract all field names of par (as saved)
    if skip_prof == 0
        skip_prof = 1; % then also skip the profiling!
        disp('Skipping profiling as you requested to use the saved set.')
        disp(' ')
    end
    disp('Calculating sample from confidence region using previously determined bounds from MAT file. Profiles are reconstructed from the MAT file as well')
elseif skip_prof == 1 && exist([filenm,'_LR.mat'],'file') == 2
    disp('Skipping profiling; using bounds from saved likelihood-region set.')
    load([filenm,'_LR']) % load the profile information from the last calc_likregion run  
    names_tmp  = fieldnames(par); % extract all field names of par (as saved)
    if length(prof_coll) ~= sum(par_sel)
        error('In the saved set, the number of fitted parameters does not equal the number of profiles,')
    end
elseif skip_prof == 2 && exist([filenm,'_LP.mat'],'file') == 2
    disp('Skipping profiling; using bounds from saved profiling set.')
    load([filenm,'_LP']) % load the profile information from the last calc_likregion run  
    names_tmp  = fieldnames(par); % extract all field names of par (as saved)
    if length(prof_coll) ~= sum(par_sel)
        error('In the saved set, the number of fitted parameters does not equal the number of profiles,')
    end    
    if ~isequal(find(pmat_tmp(:,2)==1),parnum)
        error('In the saved set, it is not the same parameters (in the same order) as required here.')
    end
else
    if skip_prof == 1
        skip_prof = 0;
        disp('There is no saved set found, so calculating profiles anyway.')
        disp(' ')
    end
    disp('Calculating profiles and sample from confidence region ... please be patient.')
end
drawnow % plot the last plot in the plot buffer before starting the analysis

pmat = packunpack(1,par,0);  % transform structure into a regular matrix
% If there are discrepancies between the par entered and the par from the
% saved set, the one from the saved set is used!

if ~isequal(pmat,pmat_tmp)
    disp(' '), warning('off','backtrace')
    if isequal(pmat(:,1:2),pmat_tmp(:,1:2)) % ah, it's just the log settings or boundaries
        warning('The log settings and/or boundaries of parameters in the saved set differs from that in the workspace at the moment. This may not hinder the analysis.')
    else
        warning('The saved parameter matrix does not (exactly) equal the one entered when calling calc_likregion. The one from the saved set is used!')
    end
    warning('on','backtrace'), disp(' ')
    fprintf('Parameter values from saved set \n');
    fprintf('=================================================================================\n');
    nfields = length(names_tmp);
    if isfield(par,'tag_fitted')
        nfields = nfields - 1; % why is this needed?
    end
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

%% Prepare for profiling

% put parameters that need to be on log scale on log10 scale (transfer needs it that way)
pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));
parshat = pmat(:,1); % this is the best-fit parameter vector
par_sel = pmat(:,2); % this is the selection vector
ind_fit = find(par_sel == 1); % index to fitted parameters
loglikmax = -1 * transfer(parshat(par_sel==1),pmat); % use transfer to obtain max log-likelihood

% For finding the bounds of the joint hypercube region, use number of
% parameters as dfs. 
chicritJ = chitable(min(5,sum(par_sel==1))); % cut-off for the 95% parameter region
chicritS = chitable(1); % cut-off for single par. CI and propagation region
if lim_out == 1
    chicritJ = chicritS + 1; % just a bit more than the inner rim
end

%% Initialise figure for plotting

[figh,ft] = make_fig(length(ind_fit),length(ind_fit)); % make figure of correct size
hold on
h1_rem = cell(length(ind_fit),1); % collect plotting handles for later use

% (re-)create the profiles on the diagonal of the multiplot
for i_p = 1:length(ind_fit) % run through the fitted parameters
    figure(figh) % make sure that the multiplot is the current plot
    h1 = subplot(length(ind_fit),length(ind_fit),(i_p-1)*length(ind_fit)+i_p); % subplot on diagonal
    hold on
    set(h1,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
    p = get(h1,'position');
    p([3 4]) = p([3 4])*1.25; % add 10 percent to width and height
    if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
        p(2) = p(2)-0.04;
    end
    set(h1, 'position', p);
    h1_rem{i_p} = h1; % remember the handle to later correct axes if needed
    
    % determine which plots get a label on which axis
    if i_p == 1
        if pmat(ind_fit(i_p),5) == 0
            ylabel(h1,[names_tmp{ind_fit(i_p)},' (log)'],ft.name,ft.label)
        else
            ylabel(h1,[names_tmp{ind_fit(i_p)}],ft.name,ft.label)
        end
        ylabel(h1,'2x loglik ratio','FontSize',12)
        set(h1,'XTickLabel',[]); % this works with old versions as well
    elseif i_p == length(ind_fit)
        if pmat(ind_fit(i_p),5) == 0
            xlabel(h1,[names_tmp{ind_fit(i_p)},' (log)'],ft.name,ft.label)
        else
            xlabel(h1,[names_tmp{ind_fit(i_p)}],ft.name,ft.label)
        end
        set(h1,'YTickLabel',[]); % remove tick labels on y-axis
    else
        set(h1,'YTickLabel',[]); % remove tick labels on y-axis
        set(h1,'XTickLabel',[]); % remove tick labels on x-axis
    end
    drawnow
end

%% Do the profiling and make plot
% The actual profiling is done by another function ... 

if skip_prof == 0 % call function that does the profiling
    opt_prof.verbose = 2; % tell calc_proflik to NOT make plots (we'll do that here later)
    % However, this setting does provide output to screen, such as CIs and
    % values for better optima if they are found.
    
    if brkprof == 2 % allow refitting if better optimum is found
        
        [par_best,XingS,prof_coll,loglik_disp] = calc_proflik(par,'all',opt_prof,opt_optim,h1_rem);
        
        if ~isempty(par_best)
            % Because we re-fitted, we need to redefine a few things!
            loglikmax = -1 * loglik_disp; % we have a new better MLL
            par  = par_best; % also needed because it will be saved in the mat file!
            pmat = packunpack(1,par,0);  % transform structure into a regular matrix
            % put parameters that need to be on log scale on log10 scale (transfer needs it that way)
            pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));
            parshat = pmat(:,1); % this is the best-fit parameter vector
        end
        
    else
        
        [par_best,XingS,prof_coll] = calc_proflik(par,'all',opt_prof,[],h1_rem);
        
        if ~isempty(par_best) % better optimum, not refitted, so have to stop!
            error('Profiling located a better optimum, no refitting set, so it makes no sense to continue.')
        end
    end
    
    load([filenm,'_LP']) % this loads the information from the profiling
    % This is needed as not all information is received as output from calc_proflik.
end

% (re-)create the profiles on the diagonal of the multiplot
for i_p = 1:length(ind_fit) % run through the fitted parameters
    figure(figh) % make sure that the multiplot is the current plot
    
    h1 = h1_rem{i_p}; % use the handle from the saved set of handles
       
    % And make a nice plot for the profiles
    plot(h1,prof_coll{i_p}(:,1),prof_coll{i_p}(:,2),'k.-')
    
    % Plot the chi-square criterion for 95% and 1 degree of freedom
    if lim_out == 0 % if we limit the space to search, chicrit is not the one for the joint CR
        plot(h1,[min(prof_coll{i_p}(:,1)) max(prof_coll{i_p}(:,1))],[chicritJ chicritJ],'k:')
    end
    plot(h1,[min(prof_coll{i_p}(:,1)) max(prof_coll{i_p}(:,1))],[chicritS chicritS],'k--')
    plot(h1,parshat(ind_fit(i_p)),0,'ko','MarkerFaceColor','y','MarkerSize',8) % plot the max lik value as a circle
    
    ylim(h1,[max(0,min(prof_coll{i_p}(:,2))) chicritJ+1]) % limit y-axis to just above the cut-off criterion
    xlim(h1,[boundscoll(i_p,1) boundscoll(i_p,2)]) % limit x-axis to saved bounds
    drawnow
    
end

%% Take an LHS (needs statistics toolbox) or normal random sample to find the joint conf. region
% Note that pmat may contain parameters that are on log-scale, and
% boundscoll too. That is correct: transfer will put it back on normal
% scale.

if nr_lhs == -1 % do NOT make a new sample, but use the saved one!
    
    disp(' ')
    disp('Using sample from saved MAT file.')
    disp(['Number of samples accepted in the total conf. region: ',num2str(size(rnd,1)),', and in inner rim: ',num2str(size(acc,1))])
    disp(' ')
    
else
    
    % Note that in pmat, log pars in first column are now on log scale!
    % Also note that sample_prof_acc contains the profile points: these are
    % included into the sample as initial points.
    nr_tot  = size(sample_prof_acc,1); % total number of parameter sets tried
    rnd     = sample_prof_acc; % start from the profiled sets, and collect the accepted sets
    n_inner = 0; % number of sets in inner rim (df=1)
    
    disp(' ')
    disp('Starting with obtaining a sample from the joint confidence region.')
    disp(['Bursts of ',num2str(burst),' samples, until at least ',num2str(nr_lhs),' samples are accepted in inner rim.'])
    
    f = waitbar(0,'Shooting for sample of joint confidence region. Please wait.','Name','calc_likregion.m');
    
    if exist('lhsdesign','file')==2 % when lhsdesign exists exists as an m-file in the path
        LHS = 1;
    else
        LHS = 0;
    end
    
    while n_inner < nr_lhs
        
        waitbar(n_inner/nr_lhs,f) % make a nice waiting bar

        if LHS == 0
            sample_lhs = rand(burst,length(ind_fit)); % uniform random sample between 0 and 1
        else
            sample_lhs = lhsdesign(burst,length(ind_fit)); % Latin-hypercube sample between 0 and 1
        end
        
        for i = 1:length(ind_fit) % go through the fitted parameters
            sample_lhs(:,i) = sample_lhs(:,i)*(boundscoll(i,2) - boundscoll(i,1))+boundscoll(i,1);
            % and change them to cover the bounds of the hypercube
        end
        
        loglik_lhs = zeros(size(sample_lhs,1),1); % pre-define for speed
        for i = 1:size(sample_lhs,1) % run through all samples in parallel
            loglik_lhs(i) = -1 * transfer(sample_lhs(i,:),pmat); % use transfer to obtain likelihood for each set
        end
        
        chi_lhs     = 2*(loglikmax - loglik_lhs); % calculate difference with the best fitting parameters that follows a chi-square
        ind_confreg = (chi_lhs<=chicritJ); % this is the index to the sets that are within the region we like to keep
        rnd         = cat(1,rnd,[sample_lhs(ind_confreg,:) chi_lhs(ind_confreg)]); % only the accepted parameter values (within the joint search region)
        % also collect the minloglik ratio as last entry in rnd
        n_inner     = n_inner + sum(chi_lhs<=chicritS); % how many are in the inner rim (df=1)
        
        % log-scale parameters are still on log-scale in this sample!
        nr_tot = nr_tot + size(sample_lhs,1);
        
    end
    
    close(f) % close the waiting bar
    
    % rnd is the accepted sample (within total 95% conf. region, df=p, or a limited region when lim_out = 1)
    rnd = sortrows(rnd,size(rnd,2)); % sort based on loglikrat (last column)
    acc = rnd(rnd(:,end)<chicritS,:); % put the sets in inner rim (df=1) into acc
    
    % save the sample in a MAT-file with the name of the script, with
    % _LR at the end. Load can be used to retrieve it.
    % also, save par (total parameter matrix), par_sel (selection
    % matrix) and bounds to make plots without redoing the fit.
    save([filenm,'_LR'],'rnd','par','par_sel','boundscoll','prof_coll','sample_prof_acc')
    % If this file already exists, this will replace it.
    
    % If we're doing a DEBtox analysis, it is a good idea to save the
    % brood-pouch delay with the sample. 
    if isfield(glo,'Tbp')
        Tbp = glo.Tbp;
        save([glo.basenm,'_LR'],'Tbp','-append')
    end
    % If we have multiple versions of the same parameter, also save that info.
    if isfield(glo,'names_sep')
        names_sep = glo.names_sep;
        save([glo.basenm,'_LR'],'names_sep','-append')
    end
    
    disp(['Total number of parameters tried (profile points have been added): ',num2str(nr_tot)])
    disp(['of which accepted in the conf. region: ',num2str(size(rnd,1)),', and in inner rim: ',num2str(size(acc,1))])
    
end

%% Add plots with the sample in all binary combinations of parameters 

if axbnds == 1 % bind axis on the bounds of the hyperbox
    minrnd = boundscoll(:,1);
    maxrnd = boundscoll(:,2);
elseif axbnds == 2 % bind axis on the bounds of the 95% region
    minrnd = min(rnd(:,1:end-1),[],1);
    maxrnd = max(rnd(:,1:end-1),[],1);
elseif axbnds == 3 % bind axis on the bounds of the inner region
    minrnd = min(acc(:,1:end-1),[],1);
    maxrnd = max(acc(:,1:end-1),[],1);
end
% adding a bit extra on the bounds seems like a good idea in all cases
minrnd = minrnd - 0.05 * (maxrnd-minrnd);
maxrnd = maxrnd + 0.05 * (maxrnd-minrnd);
    
for parnr1 = 1:length(ind_fit)-1 % go through columns
    for parnr2 = parnr1+1:length(ind_fit) % go through rows
        
        figure(figh) % make sure that the multiplot is the current plot
        h1 = subplot(length(ind_fit),length(ind_fit),(parnr2-1)*length(ind_fit)+parnr1);
        set(h1,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
        p = get(h1,'position');
        p([3 4]) = p([3 4])*1.25; % add to width and height
        if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
            p(2) = p(2)-0.04;
        end
        set(h1, 'position', p);
        hold on
        
        % Note that error ellipse is not relevant as it is not a random sample
        plot(rnd(:,parnr1),rnd(:,parnr2),'co','MarkerFaceColor','w') % plot accepted set
        plot(acc(:,parnr1),acc(:,parnr2),'ko','MarkerFaceColor','g') % plot sets within inner rim
        
        if chull == 1 && size(rnd,1)>3 % plot convex hull that approximates 95% edges of the likelihood region
            k = convhull(rnd(:,parnr1),rnd(:,parnr2));
            plot(rnd(k,parnr1),rnd(k,parnr2),'r-')
        end
        
        plot(parshat(ind_fit(parnr1)),parshat(ind_fit(parnr2)),'ko','MarkerFaceColor','y','MarkerSize',8) % plot best fit set
        axis([minrnd(parnr1) maxrnd(parnr1) minrnd(parnr2) maxrnd(parnr2)])
        
        if parnr1 == 1 % we are in the first column
            if pmat(ind_fit(parnr2),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                ylabel([names_tmp{ind_fit(parnr2)},' (log)'],ft.name,ft.label)
            else
                ylabel([names_tmp{ind_fit(parnr2)}],ft.name,ft.label)
            end
            if parnr2<length(ind_fit) % if we are not on the last row, remove axis numbers
                set(h1,'XTickLabel',[]); % this works with old versions as well
            else
                if pmat(ind_fit(parnr1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                    xlabel([names_tmp{ind_fit(parnr1)},' (log)'],ft.name,ft.label)
                else
                    xlabel([names_tmp{ind_fit(parnr1)}],ft.name,ft.label)
                end
            end
        elseif parnr2 == length(ind_fit) % we are in the last row
            if pmat(ind_fit(parnr1),5) == 0 % if parameter was fitted on log-scale, the plot is also of the log values
                xlabel([names_tmp{ind_fit(parnr1)},' (log)'],ft.name,ft.label)
            else
                xlabel([names_tmp{ind_fit(parnr1)}],ft.name,ft.label)
            end
            if parnr1>1 % if we are not on the first column, remove axis numbers
                set(h1,'YTickLabel',[]); % this works with old versions as well
            end
        else
            set(h1,'YTickLabel',[]); % this works with old versions as well
            set(h1,'XTickLabel',[]); % this works with old versions as well
        end
        
    end
end

% do something to the profile plots!
for i_p = 1:length(ind_fit) % run through the fitted parameters
    % make sure the profiles have the same x-axis as the sample
    xlim(h1_rem{i_p},[minrnd(i_p) maxrnd(i_p)])
    
    % also plot the sample in the profile plots
    hp1 = plot(h1_rem{i_p},rnd(:,i_p),rnd(:,end),'c.'); % plot accepted set
    hp2 = plot(h1_rem{i_p},acc(:,i_p),acc(:,end),'g.'); % plot sets within inner rim
    uistack(hp2,'bottom'); % move these points to background
    uistack(hp1,'bottom'); % move these points to background
end

% make an overall legend, top-right panel
h_annot = subplot(length(ind_fit),length(ind_fit),length(ind_fit)); % make an extra subplot, and remember the handle
p = get(h_annot,'position');
p([3 4]) = p([3 4])*1.2; % add to width and height
if length(ind_fit) == 2 % when plotting just 2 pars, need to move them down a bit
    p(2) = p(2)-0.04;
end
h1 = gca; % remember the current axis
set(h1,'LineWidth',1,ft.name,ft.ticks)
dim = h_annot.Position; % the position of the subplot is also the correct place for the textbox
hold on

% plot fake data points, with same symbols as in regular plots
plot(0,0,'co','MarkerFaceColor','w') % plot accepted total set
plot(0,0,'ko','MarkerFaceColor','g') % plot inner rim (df=1)
plot(0,0,'ko','MarkerFaceColor','y','MarkerSize',8)
if chull == 1 && size(rnd,1)>3 % for convex hull that approximates 95% edges of the likelihood region
    plot([0 0],[0.1 0.1],'r-')
end
plot([0 0],[0.1 0.1],'k.-') % for profile plot

% Plot the symbols for the profile cut-offs
plot([0 0],[0.1 0.1],'k:')
plot([0 0],[0.1 0.1],'k--')
plot(0,0,'ko','MarkerFaceColor','r','MarkerSize',8) % plot the max lik value as a circle

% define legend entries % removed:,['within limited set (',num2str(n_lim)' worst)']
if lim_out == 0 % if we limit the space to search, chicrit is not the one for the joint CR
    L = {['within joint CI (df=',num2str(length(ind_fit)),')'], 'within inner rim (df=1)' ...
        'best fitting set','profile likelihood',['cutoff at df=',num2str(length(ind_fit))],'cutoff at df=1'};
else
    L = {['within joint CI (limited)'], 'within inner rim (df=1)' ...
        'best fitting set','profile likelihood','within search region','cutoff at df=1'};
end
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
    h_txt = text(0.5, 1,['Likelihood-region method. Using file: ',glo.basenm, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
    set(h_txt,ft.name,ft.text); % use standard formatting for this header
end

if glo.saveplt > 0 % if we want to save the plot
    savenm = ['profile_reg_',filenm];%
    save_plot(figh,savenm);
end

disp(' ')
disp(['Time required: ' secs2hms(toc)])
