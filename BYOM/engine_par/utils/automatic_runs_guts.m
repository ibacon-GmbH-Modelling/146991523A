function [par_out,best_sel] = automatic_runs_guts(fit_tox,par,skip_sg,SEL,FEEDB,opt_optim,opt_plot,varargin)

% This is a helpful function to run a GUTS analysis in parts, based on the
% setting in fit_tox. Furthermore, it allows to run several death
% mechanisms sequentially (SD and IT, and even mixed), providing a summary
% comparison of the AICs. In general, this function will be used in
% conjunction with the parameter-space explorer (opt_optim.type = 4), but
% most of it will also work with other optimisation routines.
%
% Note: this should also work in case there are multiple data sets per
% state. However, startgrid_guts will only use the time vector of the first
% one for its derivation of search ranges for parspace.
% 
% Select what to fit with fit_tox (this is a 2-element vector). First 
% element of fit_tox is which part of the data set to use:
%   fit_tox(1) = -1  control survival (c=0) only
%   fit_tox(1) = 0   not used for GUTS
%   fit_tox(1) = 1   fit tox parameter, but when fitting, keep hb fixed; run through all
%               elements in SEL sequentially and provide a table at the end
%   fit_tox(1) = 2   fit tox parameter, but when fitting, also fit hb; run through all
%               elements in SEL sequentially and provide a table at the end
% 
% Second element of fit_tox is whether to fit or only to plot:
%   fit_tox(2) = 0   don't fit; for standard optimisations, plot results for
%               parameter values in [par], for parspace optimisations, use saved mat file.
%   fit_tox(2) = 1   fit parameters
% 
% As output:
% par_out    fitted parameter structure; when fitting tox parameters with
%            multiple configurations, the best par_out is returned
% best_sel   index for the best death mechanism
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 DATA X0mat

% Wrap the glo and glo2 structures into a single wrapper (this is an
% adaptation for the use of the parallel toolbox)
WRAP.glo  = glo;
WRAP.glo2 = glo2;

basenm_rem = glo.basenm; % remember basename as we will modify it!
best_sel   = [NaN NaN]; % this will be used to return the best death mechanism

% Remember settings so we can turn back everything at the end
glo_rem       = glo;
glo2_rem      = glo2;
X0mat_rem     = X0mat;
DATA_rem      = DATA;

opt_optim.fit = 1; % fit the parameters (1), or don't (0)
% Note: this option is placed here, rather than in the main script. Reason
% is that, perhaps counter intuitively, opt_optim.fit needs to be 1 for the
% parameter-space explorer, even when NOT fitting, otherwise the saved set
% is not loaded. The second element of fit_tox takes over as switch for
% fitting or not fitting.
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
% Note: this is set to fit parameters. Users have to set fit_tox(2) to 0 to
% use the saved set! This setting then automatically sets this
% opt_optim.ps_saved to 1 when needed.

if isempty(skip_sg) % if this variable is not input ...
    skip_sg = 0; % use startgrid for search ranges
end

id_ctrl = 0; % control is treatment that has identifier 0

switch fit_tox(1)
                
    case -1 % for fitting the survival control response only ...
        % This piece of code selects only the treatment with identifier 0
        % (zero), removes all data apart from survival, and turns fitting
        % off for all parameters but background hazard.
        
        % Some GUTS calculations can use feedbacks, with the matrix FEEDB.
        % The exact setting feedback is not relevant for the
        % control, but it must be defined to prevent errors.
        if isempty(FEEDB)
            FEEDB = [0 0 0];
        end
        if isempty(SEL) % same holds for SEL
            SEL = 1;
        end
        glo.feedb = FEEDB(1,:); % use FIRST entry in CONFIG matrix
        glo.sel   = SEL(1); % use FIRST entry in SEL matrix
        
        X0mat = X0mat(:,X0mat(1,:)==id_ctrl); % only keep controls
        DATA(:,glo.locD) = {[0]}; % remove all but the survival data from the complete data set
        % Note: this should also work in case there are multiple data sets
        % per state. However, startgrid_guts will only use the time vector
        % of the first one for its derivation of search ranges.
        
        pmat      = packunpack(1,par,[],WRAP); % turn parameter structure into a matrix
        pmat(:,2) = 0; % don't fit ANY parameters now
        par       = packunpack(2,[],pmat,WRAP); % turn parameter matrix back into a structure
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            if opt_optim.type ~= 4 % for simplex fitting
                opt_optim.fit = 0; % fit the parameters (1), or don't (0)
            else
                opt_optim.ps_saved = 1; % use saved set for parameter-space explorer (1) or not (0);
                opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
                % Fit needs to be 1 otherwise the saved set is not loaded.
            end
        else
            par.hb(2) = 1;    % but DO fit the background hazard rate
            if opt_optim.type == 4 % for the parspace explorer ...
                par = startgrid_guts(par); % replace par by estimates based on data set (as used by openGUTS)
            else % otherwise ...
                par.hb(1) = 0.01; % start from a reasonable default
            end
        end
        
        glo2.ctot = 0; % make sure that the total concentration vector to analyse only has a zero
        
        opt_plot.statsup = glo.locD; % vector with states to suppress in plotting fits (no need for damage)
        % No need for annotations or legends as we fit hb only
        opt_plot.annot   = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
        opt_plot.legsup  = 1; % set to 1 to suppress legends on fits
        glo.basenm       = [basenm_rem,'_hb']; % modify basename to make a separate PS.mat file
        glo.sel          = 1; % does not matter for fitting hb, but needs a value
        
        % optimise and plot (fitted parameters in par_out)
        par_out = calc_optim(par,opt_optim); % start the optimisation
        calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
        
        print_par(par_out,-2) % this prints out the optimised parameter values in a
        % formatted way so they can be directly copied into the main script; the
        % -2 means that only fitted parameters are displayed
        
        if opt_optim.type ~= 4 && ~isempty(varargin)
            opt_prof = varargin{1};
            calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
        end
        
    case 0 % this one is in the DEBtox version, but not in GUTS
        
        error('This option (fit_tox=0) is not available in GUTS')
         
    otherwise % for fitting tox parameters ...
        
        if fit_tox == 1 % then we fix hb to its value in par
            par.hb(2) = 0; % turn fitting for the background hazard off
        end
        par_rem = par; % remember the old par as we'll modify it in the loop below
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            if opt_optim.type ~= 4 % for simplex fitting
                opt_optim.fit = 0; % fit the parameters (1), or don't (0)
            else
                opt_optim.ps_saved = 1; % use saved set for parameter-space explorer (1) or not (0);
                opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
                % Fit needs to be 1 otherwise the saved set is not loaded.
            end
        else
            pmat = packunpack(1,par,[],WRAP); % turn output parameter structure into a matrix    
        end
        
        % Some GUTS calculations can use feedbacks, with the matrix FEEDB.
        % The exact setting feedback is not relevant for the
        % control, but it must be defined to prevent errors.
        if isempty(FEEDB) % if this variable is not in the input ...
            FEEDB = [0 0 0]; % define as no feedbacks
            feedb_chk = 0;
        else
            feedb_chk = 1;
        end
        
        % Run through settings for the death mechanism as given in SEL (SD,
        % IT and/or mixed), and feedbacks (if used) as given in FEEDB. The
        % basename, as used for naming the mat file and the output plots,
        % is adapted to identify the settings.
        MLL_coll = zeros(length(SEL),size(FEEDB,1)); % initialise matrix to catch MLLs
        par_coll = cell(length(SEL),size(FEEDB,1));  % initialise matrix to catch parameter sets
        for i = 1:length(SEL) % run through all MoAs
            for j = 1:size(FEEDB,1) % run through all feedback configuration
                
                par        = par_rem; % reset par before modifying it
                glo.sel    = SEL(i); % change global for death mechanism
                glo.feedb  = FEEDB(j,:); % change global for feedback configuration
                if feedb_chk == 0 % no feedbacks, so no need to add to basename
                    glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel)]; % create a new basenm for each death mechanism
                else
                    glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel),'_feedb',sprintf('%d',glo.feedb)];
                end
                
                switch glo.sel % make sure that right parameters are fitted
                    % This is also included in <startgrid_guts>, but would
                    % still be needed for opt_optim.type ~= 4.
                    case 1 % for SD ...
                        par.Fs(2) = 0; % never fit the threshold spread!
                    case 2 % for IT ...
                        par.bw(2) = 0; % never fit the killing rate!
                end
                % automatic creation of search ranges for GUTS tox parameters
                if fit_tox(2) == 1 && opt_optim.type == 4 && skip_sg == 0 % only when using parspace explorer for fitting
                    par = startgrid_guts(par); % calculate the search ranges and settings from data set
                end
                
                % optimise and plot (fitted parameters in par_out)
                [par_tmp,MLL] = calc_optim(par,opt_optim); % start the optimisation
                calc_and_plot(par_tmp,opt_plot); % calculate model lines and plot them
                drawnow
                MLL_coll(i,j) = MLL; % collect the MLL so we can compare them later
                par_coll{i,j} = par_tmp; % collect the par so we can output the best one
            end
        end
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            AIC_coll = nan(length(SEL),size(FEEDB,1)); % use NaNs matrix to catch MLLs
            % If you did not fit, AIC is meaningless
            if length(SEL) == 1 && size(FEEDB,1)
                best_sel = [1 1]; % indices for best configuration
                par_out  = par_coll{1,1}; % output the best par
            end
            % If a script is run to display a saved set, we can safely say
            % it's the best (that may prevent an unnecessary error).
        else
            AIC_coll  = 2*sum(pmat(:,2))+2*MLL_coll; % calculate AIC from MLL
        end
        AIC_coll  = AIC_coll - min(AIC_coll(:)); % calculate delta AIC relative to best one
        loss_coll = exp(-AIC_coll/2); % relative probability to minimise information loss
        
        % display on screen and in diary
        diary ('results.out') % collect output in the diary "results.out"
        disp(' ')
        disp('Death mech. feedbacks    MLL    delta-AIC     prob.    best')
        disp('====================================================================')
        for i = 1:length(SEL) % run through all MoAs
            for j = 1:size(FEEDB,1) % run through all feedback configuration
                str1 = sprintf('%d',SEL(i)); % string for the MoA configuration
                str2 = sprintf('%d',FEEDB(j,:)); % string for the feedbacks configuration
                fprintf(' %s           %s   %#10.2f %#10.2f %#10.2f ',str1,str2,MLL_coll(i,j),AIC_coll(i,j),loss_coll(i,j))
                if AIC_coll(i,j) == 0 % if this is the best one ...
                    fprintf('%s','     *') % give it a star
                    best_sel = [i j]; % indices for best configuration
                    par_out  = par_coll{i,j}; % output the best par
                end
                fprintf('\n')
            end
            if feedb_chk == 1
                disp('====================================================================')
            end
        end
        if feedb_chk == 0
            disp('====================================================================')
        end
        diary off  % turn off the diary function

end

% Return settings to their original value
if fit_tox(1) < 1 % when we fit tox data, it is good to keep the modified glo
    glo   = glo_rem;
end
glo2      = glo2_rem;
X0mat     = X0mat_rem;
DATA      = DATA_rem;
