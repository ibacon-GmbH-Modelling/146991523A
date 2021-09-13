function [out_conf,par] = calc_conf(par_out,opt_conf,varargin)

% Usage: [out_conf,par] = calc_conf(par_out,opt_conf,varargin)
%
% This function creates confidence intervals on model curves from the MCMC
% sample that was saved by <calc_slice.m>, or from the LHS sample that was
% saved by <calc_likregion.m>, or from the sample saved <calc_parspace.m>.
% This script only calculates the bounds; they can be plotted by a new call
% to <calc_and_plot.m> or using <plot_tktd>.
% 
% For survival analysis, we can include the sampling uncertainty into the
% confidence bounds. This is the error caused by the fact that we have a
% limited number of individuals. This is useful for model validation; e.g.,
% when using the calibrated model to predict a different exposure scenario,
% for which there are experimental data with a certain number of
% individuals. For now, this is done for Bayes only.
%
% Additionally, this function can be used to calculate sensitivities. These
% are based on the sample from a MAT file. They should be interpreted as
% the contribution of each parameter to the overall uncertainty in the
% model output.
%
% Possible options to set in a structure <opt_conf> (defined in main script
% or defaults in <prelim_checks>)
%
% Author     : Tjalling Jager 
% Date       : June 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 X0mat DATA 
global h_txt % global to allow suppressing the title text

WRAP.glo  = glo;
WRAP.glo2 = glo2;

% extract parameters from the general globals, glo and glo2
t      = glo.t;
n_X    = glo2.n_X;
filenm = glo.basenm;

% read options from structure
type_conf = opt_conf.type;     % use values from slice sampler (1), likelihood region (2), or accepted values from posterior (3) to make intervals
samerr    = opt_conf.samerr;   % include sampling error in bounds for survival data (set samerr=1)
n_samerr  = opt_conf.n_samerr; % number of sub-sampling trials for each parameter set (if samerr=1)
sens_type = opt_conf.sens;     % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
set_zero  = opt_conf.set_zero; % parameter name(s) to set to zero (usually the background hazard)
use_par_out = opt_conf.use_par_out; % set to 1 to use par_out, as entered in this function, rather than from saved set

if type_conf == 2 && samerr == 1
    samerr = 0;
    warning('off','backtrace')
    warning('Calculation of sampling error with the likelihood-region approach needs to be reconsidered, so it is turned off at the moment (still available for Bayesian analysis).')
    disp(' '), warning('on','backtrace')
end
if opt_conf.state(1) == 0
    state_sel = 1:n_X; % use all state variables for the sensitivity analysis
else
    state_sel = opt_conf.state; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
end

% The following is only needed when calculating sampling error for survival data
supp_warn = 0; % default is don't suppress warning about difference saved set at set in par_out
supp_n    = 0; % default is to derive n for sampling error from the data set for each treatment, otherwise mean is used
repls     = 1; % defauls is plot replicates
if ~isempty(varargin)
    repls = varargin{1}; % plot replicates yes/no
    if numel(varargin)>1, supp_warn = varargin{2}; end
    if numel(varargin)>2, supp_n    = varargin{3}; end
end
if ~isempty(par_out)
    pmat_tmp = packunpack(1,par_out,0,WRAP); % use this to compare to the saved one later on
else 
    pmat_tmp = -1;
end

%% Select the correct file to load, depending on the required analysis

if type_conf < 1
    error('You cannot calculate confidence intervals without specifying where to calculate them from (use opt_conf.type with a different value from zero).')
else
    [rnd,par,par_sel] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
end
if supp_warn ~= 1 && numel(pmat_tmp)>1
   par_comp(par,par_out,set_zero) % compare the saved par to par_out
end

names_tmp  = fieldnames(par); % extract all field names of par (as saved)
ind_fittag = ~strcmp(names_tmp,'tag_fitted');
names_tmp  = names_tmp(ind_fittag); % make sure that the fit tag is not in names_tmp
% Work with a copy of names, as the saved set may have different names!
% If the saved set has different names, this likely will produce an error
% elsewhere anyway. The names array is used for the legends in the
% sensitivity analysis.

loc_zero = false(length(names_tmp),1); % create logical index with zeros
if ~isempty(set_zero) % we may want to make a parameter zero (esp. background mortality)
    % allow parameters to be set to zero, such as background hazard or initial concentrations
    if ~iscell(set_zero) % for backward compatibility
        set_zero = {set_zero}; % turn it into a cell array with one element
    end
    for i = 1:length(set_zero)
        eval(['par.',set_zero{i},'(1) = 0;']); % set parameter to zero in par
        % This is likely already done in load_rnd, but it does not hurt to do it again.
        loc_zero = loc_zero == 1 | strcmp(names_tmp,set_zero{i})==1; % add this parameter to loc_zero (logical indexing)
        eval(['par_out.',set_zero{i},'(1) = 0;']); % set parameter to zero in par_out as we may want to use that one!
    end
    % rnd(:,loc_zero) = 0; % make entire random sample zero! Don't do as it may hurt the sensitivity analysis
end

%% Calculating model curves for each set in the sample

drawnow % empty plot buffer if there's something in it
n_samples = size(rnd,1); % number of sets to propagate

% initialise the matrices that will collect the output (time x scenario x samples)
X2  = cell(n_X,1); % pre-define structure
n_s = size(X0mat,2); % number of scenarios
for i = 1:n_X % run through state variables
    X2{i} = zeros(length(t),n_s,n_samples); % initialise state variable output with zeros
end

if use_par_out == 1
    par = par_out; % then we'll use the input par, rather than the saved one
    % Note: par_out must already be structured in the main script, such
    % that the fitted parameters match the ones in the saved set, etc.
end
pmat       = packunpack(1,par,0,WRAP); % transform structure from *saved set* into a regular matrix
ind_fit    = (pmat(:,2)==1); % indices to fitted parameters
ind_logfit = (pmat(:,5)==0 & pmat(:,2)==1); % indices to pars on log scale that are also fitted!
% make sure that ind_fit and ind_logfit are taken from saved set!

if size(pmat,1) ~= length(glo2.names) || sum(ind_fit) ~= size(rnd,2)
    error('There is something wrong with the parameter vector!')
end

% Start/check parallel pool
if glo2.n_cores > 0
    poolobj = gcp('nocreate'); % get info on current pool, but don't create one just yet
    if isempty(poolobj) % if there is no parallel pool ...
        parpool('local',glo2.n_cores) % create a local one with specified number of cores
    end
end

% some trickery to get parfor running ...
pmat_coll = cell(n_samples,1); % cell array to construct pmat for each sample
for k = 1:n_samples % run through all samples
    pmat(ind_fit,1) = rnd(k,:); % replace values in pmat with the k-th random sample from the MCMC
    % put parameters that need to be fitted on log scale back on normal
    % scale (as call_deri requires normal scale, in contrast to transfer.m)
    if sum(ind_logfit)>0
        pmat(ind_logfit,1) = 10.^(pmat(ind_logfit,1));
    end
    % Note: pmat is still on normal scale here, but the sample in rnd
    % contains the value on a log scale, if a parameter is fitted on log
    % scale. Call_deri needs normal scale structure par_k.
    pmat(loc_zero,1) = 0; % make parameter zero in each set of the sample!
    % this has to be done after the tranformation to normal scale!
    pmat_coll{k} = pmat; % collect it in a huge cell array
end

% some trickery to get parfor running ...
X0mat_tmp = X0mat;
glo_tmp   = glo;
Xout_coll = cell(n_s,n_samples);

parfor k = 1:n_samples % run through all samples
    pmat_tmp = pmat_coll{k}; % read it from huge matrix
    par_k = packunpack(2,0,pmat_tmp,WRAP); % transform parameter matrix into a structure
    for j = 1:n_s % run through our scenarios
        Xout_coll{j,k} = call_deri(t,par_k,X0mat_tmp(:,j),glo_tmp); % use call_deri.m to provide the output for one scenario
    end
end
clear pmat_coll % clear these large variables as they are no longer needed

for k = 1:n_samples % run through all samples
    for j = 1:n_s % run through our scenarios
        Xout = Xout_coll{j,k};
        for i = 1:n_X  % run through state variables
            X2{i}(:,j,k) = Xout(:,i); % collect state variable i into structure X2
        end
    end
end
clear Xout_coll % clear these large variables as they are no longer needed

%% Calculating intervals on the model curves

Xlo = cell(n_X,1); % pre-define structure
Xhi = cell(n_X,1); % pre-define structure

switch type_conf
    case 1 % Bayesian MCMC sample
        for i = 1:n_X % run through state variables
            Xlo{i} = prctile(X2{i},2.5,3);  % 2.5 percentile of model lines
            Xhi{i} = prctile(X2{i},97.5,3); % 97.5 percentile of model lines
        end
    case 2 % sample from likelihood region
        for i = 1:n_X % run through state variables
            Xlo{i} = min(X2{i},[],3); % minimum of model lines
            Xhi{i} = max(X2{i},[],3); % maximum of model lines
        end
    case 3 % sample from parameter space explorer
        for i = 1:n_X % run through state variables
            Xlo{i} = min(X2{i},[],3); % minimum of model lines
            Xhi{i} = max(X2{i},[],3); % maximum of model lines
        end
end

% % This is to calculate the CIs on the EC50 from the model curves
% a = X2{1}(:,1,:);
% a = squeeze(a);
% a = a ./ (ones(size(a,1),1) * a(1,:));
% EC50 = zeros(size(a,2),1);
% for i = 1:size(a,2)
%     [a_u,ind_u] = unique(a(:,i));
%     EC50(i) = interp1(a_u,t(ind_u),0.5);
% end
% switch type_conf
%     case 1
%         EC50min = prctile(EC50,2.5)
%         EC50max = prctile(EC50,97.5)
%     case 2
%         EC50min = min(EC50)
%         EC50max = max(EC50)
%     case 3
%         EC50min = min(EC50)
%         EC50max = max(EC50)
%     case 4
%         EC50min = min(EC50)
%         EC50max = max(EC50)
% end

%% Calculate and plot sensitivity of parameters over time
% Here, I use the sample that is also used for the confidence bounds on the
% model curves. Note that the correlations between parameters will also
% affects the sensitivities calculated from them.

if sens_type>0 % only if we ask for it
    for ii = 1:length(state_sel) % run through the states that we like to see
        
        ind_rnd = find(pmat(:,2)==1); % index in pmat for the parameters in the sample rnd
        corr_coll = nan(length(t),length(ind_rnd),n_s); % initialise correlations matrix with NaNs
        
        for k = 1:n_s % run through scenarios
            for i = 1:length(t) % run through time points
                for j = 1:length(ind_rnd) % run through parameters in the sample
                    
                    % it is good to check if there is meaningful variation
                    % in the calculated state variable; otherwise, the
                    % corrcoef will just reflect the accuracy of the ODE
                    % solver
                    Xi = X2{state_sel(ii)}(i,k,:); % these are the values for the state/time/scenario that we will correlate
                    Xi_cv = std(Xi)/mean(Xi);
                    if isnan(Xi_cv) || Xi_cv < 1e-6
                        tmp = [0 0];
                    else
                        switch sens_type
                            case 1
                                tmp = corrcoef(Xi,rnd(:,j)); % corr. parameter value and state
                            case 2
                                tmp = corrcoef(Xi./X2{state_sel(ii)}(i,1,:),rnd(:,j));% corr. parameter value and effect relative to treatment 1
                            case 3
                                if i < length(t) % we look at change, so ignore last time point
                                    rel_state = (X2{state_sel(ii)}(i+1,k,:)-Xi)./Xi;
                                    tmp = corrcoef(rel_state,rnd(:,j)); % corr. parameter value and effect relative to state
                                else
                                    tmp = [NaN NaN; NaN NaN];
                                end
                        end
                    end
                    corr_coll(i,j,k) = tmp(2); % remember the correlation that matters
                end
            end
        end
        
        % make one plot with subplots for each scenario (i.e., concentration)
        n = ceil(sqrt(n_s));
        m = ceil(n_s/n);
        plotcol = 'krbgcm'; % 6 colors should be enough ...
        
        figh = make_fig(m,n); % create figure window of correct size
        
        for k = 1:n_s % run through scenarios
            
            if n_s > 1
                g = subplot(m,n,k);
            end
            set(gca,'LineWidth',1,'FontSize',12) % adapt axis formatting
            if m > 1 % for 1 row, making the plot larger is no fun
                p = get(g,'position');
                p([3 4]) = p([3 4])*1.10; % add 10 percent to width and height
                set(g, 'position', p);
                title([glo.leglab1,num2str(X0mat(1,k)),' ',glo.leglab2],'FontSize',10)
            end
            
            hold on
            
            for j = 1:length(ind_rnd) % run through parameters in the sample
                if j <= 6
                    plot(t,corr_coll(:,j,k),[plotcol(j),'-'],'LineWidth',2)
                else % for more than 6 parameters, use dotted line
                    plot(t,corr_coll(:,j,k),[plotcol(j-6),':'],'LineWidth',2)
                end
            end
            
            ax = gca; % handle to current axis system
            if k>(n*(m-1)) % only put xlabel on bottom row
                xlabel(glo.xlab,'FontSize',12)
            else
                % ax.XTickLabel = {[]}; % this is okay for new Matlab versions
                set(ax,'XTickLabel',[]); % this works with old versions as well
            end
            if k==1 || (k-1)/n == floor((k-1)/n) % only put y labels on first column
                switch sens_type
                    case 1
                        ylabel(['correlation to state ',num2str(state_sel(ii))],'FontSize',12)
                    case 2
                        ylabel(['corr. to state/control ',num2str(state_sel(ii))],'FontSize',12)
                    case 3
                        ylabel(['corr. rel. change of state ',num2str(state_sel(ii))],'FontSize',12)
                end
            else
                % ax.YTickLabel = {[]};
                set(ax,'YTickLabel',[]);
            end
            
            plot([t(1) t(end)],[0 0 ],'k:','LineWidth',2)
            axis([t(1) t(end) -1 +1]) % set all axis within the multiplot the same
            
        end
        
        legend(names_tmp{ind_rnd})
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
        h_txt = text(0.5, 1,['Contribution to uncertainty state ',num2str(state_sel(ii))],'HorizontalAlignment','center','VerticalAlignment', 'top');
        
        if glo.saveplt > 0 % if we want to save the plot
            savenm = ['confsense_',filenm,'_state_',num2str(ii)];%
            save_plot(figh,savenm,h_txt);
        end
    end
end

%% Include sampling uncertainty in intervals for survival data
% The script <calc_and_plot> will plot the intervals including the sampling
% uncertainty as thinner broken lines. For MCMC samples, I just calculate a
% number of multinomial trials (<n_samerr>) from each element of the
% sample. From the entire forest of curves (<n_samerr> times number of
% elements in the sample), percentiles are taken.
% 
% For likelihood-region sample, the derivation is more awkward, owing to
% the awkward definition of CRs and CIs. I do the same thing, but this is
% questionable (though it yields comparable results to Bayes). I have tried
% a different way, but this overestimates the prediction interval: for each
% element of the sample, <n_samerr> multinomial trials will be run. First,
% percentiles are taken on these trials, and these will be collected for
% all elements of the sample. The final intervals will then be the
% minumum-maximum of these values. At some point, this needs to be replaced
% by the 'validation profile likelihood' of Kreutz et al.

XloS = -1;
XhiS = -1;
XloD = -1;
XhiD = -1;

if samerr == 1 && isfield(glo,'locS') && numel(DATA{1,glo.locS})>1 && exist('mnrnd','file')==2 % only if mnrnd is in the path (statistics toolbox needed)
    
    % initialise the matrices that will collect the output (time x scenario x samples)
    XloS = nan(length(t),n_s); % initialise with zeros
    XhiS = nan(length(t),n_s); % initialise with zeros
    % initialise the matrices that will collect the output for deaths (time x scenario x samples)
    XloD = nan(length(t),n_s); % initialise with zeros
    XhiD = nan(length(t),n_s); % initialise with zeros
    
    % Need to catch cases where the number of individuals in the test is
    % not know. For example, when using calc_doseresp (and many
    % concentrations are simulated that are not in the data set), or when
    % simulations are made for extra scenarios.
    cd1   = DATA{1,glo.locS}(1,2:end); % all treatments in the data set
    cd2   = unique(cd1); % unique treatments in the data set
    n_rep = floor(mean(DATA{1,glo.locS}(2,2:end))); % mean number of individuals in each replicate
    n_mn  = sum(DATA{1,glo.locS}(2,2:end))/length(cd2); % mean number of individuals per treatment    
    
    % this is for the progress displayer
    disp('Calculating bounds including sampling error, please be patient ...')
    fprintf('Working on scenario (of %d): ',n_s)
    
    n_check = 0; % check whether n=0 occurred
    for j = 1:n_s % run through our scenarios
        
        fprintf('%d ',j)
        if j/20 == floor(j/20), fprintf('\n'), end % create a new line once in a while
        
        % first, derive number of individuals from the data set
        % this works as long as there is 1 data set with survival
        if supp_n == 1 % for dose-response plots, I always use the mean over all treatments/replicates
            % otherwise, you get nasty error bounds, depending on whether
            % you calculate a concentration that is in the data set or not
            if repls == 1, n = n_rep; 
            else n = n_mn; end
        else
            if repls == 1 % when plotting individual replicates ...
                n = floor(mean(DATA{1,glo.locS}(2,find(DATA{1,glo.locS}(1,:)==X0mat(1,j)))));
                % NOTE: I use mean here, so if there are multiple entries with the
                % same exposure concentration, the sampling error represents
                % the individual replicates with mean nr. individuals.
                if n == 0 || isnan(n), n = n_rep; n_check = 1; end % if there is no answer, use mean
            else % if plotting the mean, the sample error is based on the sum of the replicates
                n = sum(DATA{1,glo.locS}(2,find(DATA{1,glo.locS}(1,:)==X0mat(1,j))));
                if n == 0 || isnan(n), n = n_mn; n_check = 1; end % if there is no answer, use mean
            end
        end
        
        Scoll = []; % initialise for collection of samplings for individuals
        Dcoll = []; % initialise for collection of samplings for individuals
        
        for k = 1:n_samples % run through all samples
            
            S = X2{glo.locS}(:,j,k); % survival probabilities
            p = [-diff(S);1+sum(diff(S))]; % multinomial probabilities
            if any(p) < 0 % sometimes one element gets slightly negative ...
                p = max(0,p); % set them to zero ... and ...
                p = p./sum(p); % normalise again to make sure the sum is one
            end
            
            Scollk = nan(length(t),n_samerr);
            Dcollk = nan(length(p),n_samerr);
            for i = 1:n_samerr % random samples for survivors
                deaths = mnrnd(n,p); % take a sample from multinomial distribution
                Ssim = [n n-cumsum(deaths)]/n; % and translate to survival prob.
                Ssim(end) = []; % last interval is after end of test, so remove
                Scollk(:,i) = Ssim';
                Dcollk(:,i) = deaths'; % for the death plots!
            end
            
            % Modified the calculation of bounds; commenting these parts
            % out seems to lead to more realistic results for the likregion.
%             if type_conf == 1 % collect all results
                Scoll = cat(2,Scoll,Scollk); % add the results for each item in the sample to the list
                Dcoll = cat(2,Dcoll,Dcollk); % add the results for each item in the sample to the list
%             else % take percentiles for the sampling errors, so that we can take min-max later on!
%                 Scoll = cat(2,Scoll,[prctile(Scollk,2.5,2) prctile(Scollk,97.5,2)]); % add the percentiles for each item in the sample to the list
%                 Dcoll = cat(2,Dcoll,[prctile(Dcollk,2.5,2) prctile(Dcollk,97.5,2)]); % add the percentiles for each item in the sample to the list
%             end
             
        end
        
%         if type_conf == 1
            XloS(:,j) = prctile(Scoll,2.5,2);  % 2.5 percentile of model lines
            XhiS(:,j) = prctile(Scoll,97.5,2); % 97.5 percentile of model lines
            XloD(:,j) = prctile(Dcoll,2.5,2);  % 2.5 percentile of model lines
            XhiD(:,j) = prctile(Dcoll,97.5,2); % 97.5 percentile of model lines
%         else
%             XloS(:,j) = min(Scoll,[],2); % minimum of model lines
%             XhiS(:,j) = max(Scoll,[],2); % maximum of model lines
%             XloD(:,j) = min(Dcoll,[],2); % minimum of model lines
%             XhiD(:,j) = max(Dcoll,[],2); % maximum of model lines
%         end
        
    end
    
    fprintf('\n') % end the line opened by the progress monitor
    
    if n_check == 1
        disp(' '), warning('off','backtrace')
        warning('For some calculated scenarios, there was no information on number of individuals. A mean of all treatments was applied.')
        warning('off','backtrace'), disp(' ')
    end

end

disp(['Time required: ' secs2hms(toc)])

% collect all output into a structure
out_conf    = cell(8,1);
out_conf{1} = Xlo;
out_conf{2} = Xhi;
out_conf{3} = XloS;
out_conf{4} = XhiS;
out_conf{5} = XloD;
out_conf{6} = XhiD;
out_conf{7} = X2;
out_conf{8} = type_conf;

