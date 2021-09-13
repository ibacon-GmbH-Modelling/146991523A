function [par_out,best_MoaFb] = automatic_runs_debtox2019(fit_tox,par,ind_tox,skip_sg,MOA,FEEDB,opt_optim,opt_plot,varargin)

% This is a helpful function to run a DEBtox analysis in parts, based on
% the setting in fit_tox, and to sequentially run through a series of model
% configurations. In general, this function will be used in conjunction
% with the parameter-space explorer (opt_optim.type = 4), at least for
% fitting the toxicity data, but most of it will also work with other
% optimisation routines.
% 
% Select what to fit with fit_tox (this is a 2-element vector). First 
% element of fit_tox is which part of the data set to use:
%   fit_tox(1) = -2  analyse control data (compare segular/solvent controls
%                    and controls across data sets)
%   fit_tox(1) = -1  fit/show control survival only
%   fit_tox(1) = 0   fit/show controls for growth/repro only, but not for survival
%   fit_tox(1) = 1   fit/show all treatments, but, when fitting, keep all control parameters fixed; 
%               run through all elements in MOA and FEEDB sequentially and 
%               provide a table at the end (plots are made incl. control)
% 
% Second element of fit_tox is whether to fit or only to plot:
%   fit_tox(2) = 0   don't fit; for standard optimisations, plot results for
%               parameter values in [par], for parspace optimisations, use saved mat file.
%   fit_tox(2) = 1   fit parameters
% 
% Third element is what to use as control (fitted for fit_tox(1) = -1 or 0)
% (if this element is not present, only regular control will be used)
%   fit_tox(3) = 1   use regular control only (identifier 0)
%   fit_tox(3) = 2   use solvent control only (identifier 0.1)
%   fit_tox(3) = 3   use both regular and control (identifier 0 and 0.1)
% 
% As output:
% par_out    fitted parameter structure; when fitting tox parameters with
%            multiple configurations, the best par_out is returned
% best_MoaFb indices for the best pMoA and feedback configuration
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

% IF THERE ARE TWO DATA SETS, MAKE SURE THAT SECOND DATA SET HAS
% IDENTIFIERS STARTING WITH 100 FOR CONTROL AND 100.1 FOR THE SOLVENT
% CONTROL!!

global glo glo2 DATA W X0mat

% Wrap the glo and glo2 structures into a single wrapper (this is an
% adaptation for the use of the parallel toolbox)
WRAP.glo  = glo;
WRAP.glo2 = glo2;

basenm_rem = glo.basenm; % remember basename as we will modify it!
best_MoaFb = [NaN NaN]; % this will be used to return the best MoA-feedback configuration

% Remember settings so we can turn back everything at the end
glo_rem       = glo;
glo2_rem      = glo2;
X0mat_rem     = X0mat;
DATA_rem      = DATA;
W_rem         = W;
par_rem       = par;

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

if isempty(skip_sg) % if this variable is not in the input
    skip_sg = 0; % use startgrid for search ranges
end
if isempty(MOA) % if this variable is not in the input
    MOA   = [1 0 0 0 0];
end
if length(MOA) == 4
    MOA(5) = 0; % make sure the length of MOA is 5
end
if isempty(FEEDB) % if this variable is not in the input
    FEEDB = [0 0 0 0];
end
% MOA and FEED need be defined for fitting the controls, but settings don't matter

id_solvent = 0.1; % identifier for solvent control (control should be zero)
% this could become a global or part of an option structure in future versions

if length(fit_tox) < 3 % it might be missing ...
    fit_tox(3) = 1; % by default, only use regular control (identifier 0)
end
% Decide what to use as controls
switch fit_tox(3)
    case 1 % only use regular control
        id_ctrl = [0];
    case 2 % only use solvent control
        id_ctrl = [id_solvent];
    case 3 % use both regular and solvent control
        id_ctrl = [0 id_solvent];
end

switch fit_tox(1)
    
    case -2 % compare controls (or solvent controls) across data sets (likelihood-ratio test)
        
        % Never show iterations for these control fits
        opt_optim.it = 0; % show iterations of the simplex optimisation (1, default) or not (0)
        
        % TEMP TEMP this may be better in the main script or prelim_checks???
        if isempty(glo.names_sep) % if empty, use a default for this analysis
            names_sep = {'f';'L0'}; % names of parameters that can differ between data sets (for all but f: only when fit mark is 1)
        else 
            names_sep = glo.names_sep; % use names for separate parameters as given in glo
        end
        glo.names_sep = {}; % but always start with empty global for separate name; we fit individual data sets first
        
        % The exact setting for moa and feedb is not relevant for the
        % control, but they must be defined to prevent errors.
        glo.moa   = MOA(1,:);   % use FIRST entry in MOA matrix
        glo.feedb = FEEDB(1,:); % use FIRST entry in FEEDB matrix
        
        ns       = size(DATA_rem,2); % number of states
        nd       = size(DATA_rem,1); % number of data sets
        
        flds_rem = 0; % counter for number of removed fields
        if nd > 1 % if there are more data sets
            % Make sure that any additional parameters for other data sets
            % are NOT used initially! They will be generated again when
            % testing for common parameters across data sets.
            if par.f(2) == 1
                warning('You are fitting <f> in the *first* data set ... this could easily lead to nonsense ...')
            end
            for i_set = 1:nd-1 % run through data sets (not first one)
                for i_sep = 1:length(names_sep) % run through extra parameter names for separate sets
                    str_tst = [names_sep{i_sep},num2str(i_set)];
                    if isfield(par,str_tst)
                        par      = rmfield(par,str_tst); % remove field!
                        flds_rem = flds_rem + 1; % count number of removed fields
                    end
                end
            end
            names      = fieldnames(par); % extract all field names of par (global)
            glo2.names = names; % make sure glo2 has the correct pruned names set
        end
        if flds_rem > 0
            ind_tox = ind_tox - flds_rem;
            warning('Removed all already-defined extra parameters for additional data sets from par. Check results carefully! I assume that they were added BEFORE the definition of ind_tox.')
        end
        
        DATA_ctrl = cell(nd,ns);      % create an empty data set with single data set
        W_ctrl    = cell(nd,ns);      % create an empty weight set with single data set
        nosurv    = 0; % flag for if we don't have survival data
        % Note: at this moment, if there are no survival data for 1 set,
        % survival is not used at all.
        
        % create a new data set with only control data, for all sets
        for i_set = 1:nd % run through data sets
            for i = 1:ns % run through all states
                if ~isempty(DATA_rem{i_set,i}) && numel(DATA_rem{i_set,i}) > 1 % if there is data ...
                    id = DATA_rem{i_set,i}(1,2:end); % identifiers
                    
                    % find indices for controls, also when there is a
                    % regular AND a solvent control
                    [~,ind_id_c]  = ismember(id,id_ctrl+(i_set-1)*100); % find where controls are in id
                    ind_id_c = find(ind_id_c>0);
                    
                    DATA_ctrl{i_set,i} = [DATA_rem{i_set,i}(:,1),DATA_rem{i_set,i}(:,ind_id_c+1)];
                    W_ctrl{i_set,i} = W_rem{i_set,i}(:,ind_id_c);
                    
                else
                    DATA_ctrl{i_set,i} = 0;
                    if i == glo.locS % then survival is one of the empty states
                        nosurv = 1; % set flag
                    end
                end
            end
        end
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            if opt_optim.type ~= 4 % for simplex fitting
                opt_optim.fit = 0; % fit the parameters (1), or don't (0)
            else
                opt_optim.ps_saved = 1; % use saved set for parameter-space explorer (1) or not (0);
                opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
                % Fit needs to be 1 otherwise the saved set is not loaded.
            end
        end
        
        opt_plot.statsup = [glo.locD]; % vector with states to suppress in plotting fits
        if nosurv == 1 % then we don't have survival data
            par.hb(2) = 0; % don't fit the background hazard rate
            opt_plot.statsup = [glo.locD glo.locS]; % vector with states to suppress in plotting fits
        end
        WRAP.glo2 = glo2; % make sure WRAP is updated as glo2 might have changed
        pmat   = packunpack(1,par,[],WRAP); % turn parameter structure into a matrix
        pmat(ind_tox:end,2) = 0; % don't fit the toxicity parameters now
        par    = packunpack(2,[],pmat,WRAP); % turn parameter matrix back into a structure
        nr_fit = sum(pmat(:,2));
        
        % No need for annotations (and perhaps legends) as we fit basic parameters only
        opt_plot.annot  = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
        opt_plot.legsup = 0; % set to 1 to suppress legends on fits
        
        % =================================================================
        % First, fit all controls separately
        % =================================================================
        par_coll   = nan(length(id_ctrl),nr_fit,nd); % matrix to collect parameter fits
        mll_coll   = nan(length(id_ctrl),nd); % matrix to collect MLL for each fit
        X0mat_coll_c = []; % collect X0mat for all regular controls
        X0mat_coll_s = []; % collect X0mat for all solvent controls
        
        for i_set = 1:nd % run through data sets (if there are multiple)
            for i_c = 1:length(id_ctrl) % run through regular/solvent control (if they're there)
                
                DATA     = DATA_ctrl(i_set,:); % only use data for this set
                W        = W_ctrl(i_set,:); % only use weights for this set
                glo2.n_D = 1; % tell everyone there is only one data set now!
                
                [~,loc_i]  = ismember(id_ctrl(i_c)+(i_set-1)*100,X0mat_rem(1,:)); % find where this control is in X0mat
                X0mat      = X0mat_rem(:,loc_i); % use only that column for new global X0mat
                
                if id_ctrl(i_c) == 0
                    X0mat_coll_c = [X0mat_coll_c,X0mat]; % add this temporary X0mat to the collection
                else
                    X0mat_coll_s = [X0mat_coll_s,X0mat]; % add this temporary X0mat to the collection
                end
                
                glo2.ctot  = X0mat(1,:); % make sure that the total concentration vector to analyse only has the identifier for the control
                glo.basenm = [basenm_rem,'_ctrl_set_',num2str(i_set)]; % modify basename to make a separate PS.mat file
                [par_out,FVAL] = calc_optim(par,opt_optim); % start the optimisation (don't plot)
                % calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
                
                WRAP.glo2 = glo2; % make sure WRAP is updated as glo2 might have changed
                pmat_tmp = packunpack(1,par_out,[],WRAP); % turn parameter structure into a matrix
                par_coll(i_c,:,i_set) = pmat_tmp(pmat_tmp(:,2)==1,1); % collect fitted parameters
                mll_coll(i_c,i_set)   = FVAL; % collect MLL
                
                if id_ctrl(i_c) == 0
                    disp(['Control fit data set ',num2str(i_set),', regular control'])
                else
                    disp(['Control fit data set ',num2str(i_set),', solvent control'])
                end
                print_par(par_out,-2) % this prints out the optimised parameter values in a
                % formatted way so they can be directly copied into the main script; the
                % -2 means that only fitted parameters are displayed
                disp(' ')
            end
        end
        
        % Chi2 as a table so we don't need statistics toolbox!
        crit_table = [3.8415
            5.9915
            7.8147
            9.4877
            11.07
            12.592
            14.0671]; % critical values (0.05) from chi2 at different DFs
        
        % =================================================================
        % Next, fit control and solvent control together for each data set
        % =================================================================
        par_coll_cs = nan(nr_fit,nd); % collect parameters from fits
        mll_coll_cs = nan(nd,1); % collect MLL from fits
        if fit_tox(3) == 3 % only if both controls are requested
            for i_set = 1:nd % run through data sets
                DATA = {}; % clear previous data set
                W    = {}; % clear previous data set
                
                for i = 1:ns % run through all states
                    
                    id = DATA_ctrl{i_set,i}(1,2:end); % identifiers
                    [~,ind_id_c]  = ismember(id,id_ctrl(1)+(i_set-1)*100); % find where regular controls are in id
                    ind_id_c = find(ind_id_c>0);
                    [~,ind_id_s]  = ismember(id,id_ctrl(2)+(i_set-1)*100); % find where solvent controls are in id
                    ind_id_s = find(ind_id_s>0);
                    
                    if ~isempty(DATA_ctrl{i_set,i}) && numel(DATA_ctrl{i_set,i}) > 1 % if there is data ...
                        DATA{1,i} = [DATA_ctrl{i_set,i}(:,1),DATA_ctrl{i_set,i}(:,ind_id_c+1)];
                        DATA{2,i} = [DATA_ctrl{i_set,i}(:,1),DATA_ctrl{i_set,i}(:,ind_id_s+1)];
                        W{1,i} = W_ctrl{i_set,i}(:,ind_id_c);
                        W{2,i} = W_ctrl{i_set,i}(:,ind_id_s);
                    else
                        DATA{1,i} = 0;
                        DATA{2,i} = 0;
                        if i == glo.locS % then survival is one of the empty states
                            nosurv = 1; % set flag
                        end
                    end
                end
                
                glo2.n_D   = 2; % tell everyone there are TWO data sets now!
                [~,loc_i]  = ismember(id_ctrl+(i_set-1)*100,X0mat_rem(1,:)); % find where this control is in X0mat
                X0mat      = X0mat_rem(:,loc_i); % use only those columns
                
                glo2.ctot  = X0mat(1,:); % make sure that the total concentration vector to analyse only has the identifier for the control
                glo.basenm = [basenm_rem,'_allctrl_set_',num2str(i_set)]; % modify basename to make a separate PS.mat file
                [par_out,FVAL] = calc_optim(par,opt_optim); % start the optimisation (don't plot)
                
                if nd == 1 % if there is only 1 data set, plot this fit
                    calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
                end
                
                WRAP.glo2 = glo2; % make sure WRAP is updated as glo2 might have changed
                pmat_tmp   = packunpack(1,par_out,[],WRAP); % turn parameter structure into a matrix
                par_coll_cs(:,i_set) = pmat_tmp(pmat_tmp(:,2)==1,1); % collect fitted parameters
                mll_coll_cs(i_set) = FVAL;
                
                disp(['Combined controls (regular and solvent) fit data set ',num2str(i_set)])
                print_par(par_out,-2) % this prints out the optimised parameter values in a
                % formatted way so they can be directly copied into the main script; the
                % -2 means that only fitted parameters are displayed
                disp(' ')
                
            end
        end
        
        DATA     = DATA_ctrl; % use all data for controls again
        W        = W_ctrl; % use all weights for controls again
        glo2.n_D = nd; % tell everyone there is the complete number of data sets again
        
        % =================================================================
        % Fit all data sets simultaneously, with the same set of parameters
        % or with selected separate parameters; this part fits all regular
        % controls together, and next all solvent controls together.
        % =================================================================
        mll_coll_d = nan(length(id_ctrl),2);
        if nd > 1 % if there are more data sets
            
            % Create a second parameter set, with extra parameter entries
            % for parameters that are allowed to differ between data sets
            par2 = par; % copy par
            for i_set = 1:nd-1 % run through data sets (not first one)
                for i_sep = 1:length(names_sep) % run through extra parameter names for separate sets
                    if isfield(par,names_sep{i_sep})
                        eval(['par2.',[names_sep{i_sep},num2str(i_set)],' = par.',names_sep{i_sep},';']); % create extra parameter, as copy of first one
                        if strcmp('f',names_sep{i_sep})
                            eval(['par2.f',num2str(i_set),'(2) = 1;']); % fit it! (needed as first one is usually fixed to 1)
                        end
                    end
                end
            end
            names2     = fieldnames(par2); % extract all field names of par (global)
            % ind_fittag = ~strcmp(names2,'tag_fitted');
            % names2     = names2(ind_fittag); % make sure that the fit tag is not in names
            names      = glo2.names; % read original names from global
            
            glo2.names = names2; % put the modified one in the global
            WRAP.glo2 = glo2; % make sure WRAP is updated as glo2 might have changed
            pmat2      = packunpack(1,par2,[],WRAP); % turn parameter structure into a matrix
            nr_fit2    = sum(pmat2(:,2)); % extract the number of fitted parameters
            glo2.names = names; % put the original back into the global
            
            for i_c = 1:length(id_ctrl) % run through regular/solvent control
                
                if id_ctrl(i_c) == 0
                    X0mat      = X0mat_coll_c;
                    glo.basenm = [basenm_rem,'_regular_ctrl_set_all']; % modify basename to make a separate PS.mat file
                else
                    X0mat      = X0mat_coll_s;
                    glo.basenm = [basenm_rem,'_solvent_ctrl_set_all']; % modify basename to make a separate PS.mat file
                end
                glo2.ctot = X0mat(1,:); % make sure that the total concentration vector to analyse only has a zero
                
                glo2.names = names; % make sure the global has the original names
                [par_out,FVAL] = calc_optim(par,opt_optim); % start the optimisation (don't plot)
                % calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
                mll_coll_d(i_c,1) = FVAL;
                
                WRAP.glo2 = glo2; % make sure WRAP is updated as glo2 might have changed
                pmat_tmp   = packunpack(1,par_out,[],WRAP); % turn parameter structure into a matrix
                %         par_coll(:,i_set) = pmat_tmp(pmat_tmp(:,2)==1,1); % collect fitted parameters
                %         mll_coll(i_set) = FVAL;
                
                if id_ctrl(i_c) == 0
                    disp(['Control fit across all data sets simultaneously, regular control '])
                else
                    disp(['Control fit across all data sets simultaneously, solvent control '])
                end
                print_par(par_out,-2) % this prints out the optimised parameter values in a
                % formatted way so they can be directly copied into the main script; the
                % -2 means that only fitted parameters are displayed
                
                disp(' ')
                disp('Compare to parameter values from individual fit on different sets:')
                fit_names = glo2.names(pmat_tmp(:,2)==1);
                fprintf('%-5s ','Set')
                for j = 1:nd
                    fprintf('%10.0f',j)
                end
                fprintf('\n')
                disp('-------------------------------------------------------')
                for i = 1:nr_fit % run through fitted parameters
                    fprintf('%-5s ',fit_names{i})
                    for j = 1:nd
                        fprintf('%#10.4g',par_coll(i_c,i,j))
                    end
                    fprintf('\n')
                end
                disp('-------------------------------------------------------')
                
                if id_ctrl(i_c) == 0
                    glo.basenm = [basenm_rem,'_regular_ctrl_set_all_seppar']; % modify basename to make a separate PS.mat file
                else
                    glo.basenm = [basenm_rem,'_solvent_ctrl_set_all_seppar']; % modify basename to make a separate PS.mat file
                end
                
                glo2.names = names2; % put the modified names in the global
                glo.names_sep = names_sep; % and notify call_deri that we have separate parameters to take care of
                [par_out,FVAL] = calc_optim(par2,opt_optim); % start the optimisation (don't plot)
                % calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
                mll_coll_d(i_c,2) = FVAL;
                
                if id_ctrl(i_c) == 0
                    disp(['Control fit all data sets simultaneously, separate parameters, regular control '])
                else
                    disp(['Control fit all data sets simultaneously, separate parameters, solvent control '])
                end
                
                fprintf('Separately fitted per data set:')
                for i_sep = 1:length(names_sep) % run through extra parameter names
                    fprintf(' %s',names_sep{i_sep})
                    if i_sep < length(names_sep)
                        fprintf(',')
                    end
                end
                fprintf('\n')
                
                print_par(par_out,-2) % this prints out the optimised parameter values in a
                % formatted way so they can be directly copied into the main script; the
                % -2 means that only fitted parameters are displayed
                disp(' ')
                
                glo.names_sep = {}; % set it back to empty again
                glo2.names    = names; % make sure the global has the original names
                
            end
            
            glo2.names = names; % put it in the global
            
            % =============================================================
            % Fit regular+solvent together, for each data set separately.
            % =============================================================
            if fit_tox(3) == 3 % if there are multiple controls ...
                % we also can compare individual fits on all data sets
                % separately to a simultaneous fit (so lumping regular and
                % solvent control). Note that this is different from the test
                % for lumping solvent and regular control for each data set: in
                % that test, separate data sets were made for the solvent and
                % regular controls, while here, they are kept within their own
                % set. We are not comparing solvent to regular, but lumping
                % them to compare individual fits on each data set to
                % simultaneous fits.
                
                par_coll_cs2 = nan(nr_fit,nd); % matrix to collect parameter fits
                mll_coll_cs2 = nan(nd,1); % matrix to collect MLL for fits
                
                for i_set = 1:nd % run through data sets
                    
                    DATA     = DATA_ctrl(i_set,:); % only use data for this set
                    W        = W_ctrl(i_set,:); % only use weights for this set
                    glo2.n_D = 1; % tell everyone there is only one data set now!
                    
                    [~,loc_i]  = ismember(id_ctrl+(i_set-1)*100,X0mat_rem(1,:)); % find where the controls are in X0mat
                    X0mat      = X0mat_rem(:,loc_i); % use only those columns
                    
                    glo2.ctot  = X0mat(1,:); % make sure that the total concentration vector to analyse only has the identifier for the control
                    glo.basenm = [basenm_rem,'_allctrl2_set_',num2str(i_set)]; % modify basename to make a separate PS.mat file
                    [par_out,FVAL] = calc_optim(par,opt_optim); % start the optimisation (don't plot)
                    % calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
                    
                    WRAP.glo2 = glo2; % make sure WRAP is updated as glo2 might have changed
                    pmat_tmp = packunpack(1,par_out,[],WRAP); % turn parameter structure into a matrix
                    par_coll_cs2(:,i_set) = pmat_tmp(pmat_tmp(:,2)==1,1); % collect fitted parameters
                    mll_coll_cs2(i_set)   = FVAL;
                    
                    disp(['Control fit (regular and solvent in same data set, so same residual s.e.) data set ',num2str(i_set)])
                    print_par(par_out,-2) % this prints out the optimised parameter values in a
                    % formatted way so they can be directly copied into the main script; the
                    % -2 means that only fitted parameters are displayed
                    disp(' ')
                end
                
            end
            
            DATA     = DATA_ctrl; % use all data for controls again
            W        = W_ctrl; % use all weights for controls again
            glo2.n_D = nd; % tell everyone there are more data sets again
            
            % =============================================================
            % Next, fit ALL data together: regular+solvent, across all data
            % sets.
            % =============================================================
            if fit_tox(3) == 3 % only when there are two controls to lump
                
                X0mat      = [X0mat_coll_c X0mat_coll_s]; % combine controls
                glo.basenm = [basenm_rem,'_regular_allctrl_set_all']; % modify basename to make a separate PS.mat file
                glo2.ctot  = X0mat(1,:); % make sure that the total concentration vector to analyse only has a zero
                
                [par_out,FVAL] = calc_optim(par,opt_optim); % start the optimisation (don't plot)
                % calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
                mll_coll_d_all = FVAL;
                disp(['Control fit (regular and solvent controls lumped) across all data sets simultaneously'])
                
                print_par(par_out,-2) % this prints out the optimised parameter values in a
                % formatted way so they can be directly copied into the main script; the
                % -2 means that only fitted parameters are displayed
                disp(' ')
                
                glo.basenm = [basenm_rem,'_regular_allctrl_set_all_seppar']; % modify basename to make a separate PS.mat file
                glo2.names = names2; % put the modified names in the global
                glo.names_sep = names_sep; % and notify call_deri that we have separate parameters to take care of
                [par_out,FVAL] = calc_optim(par2,opt_optim); % start the optimisation (don't plot)
                
                calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
                
                mll_coll_d_sep = FVAL;
                
                disp(['Control fit (regular and solvent controls lumped) across all data sets simultaneously, separate parameters'])
                
                fprintf('Separately fitted per data set:')
                for i_sep = 1:length(names_sep) % run through extra parameter names
                    fprintf(' %s',names_sep{i_sep})
                    if i_sep < length(names_sep)
                        fprintf(',')
                    end
                end
                fprintf('\n')
                print_par(par_out,-2) % this prints out the optimised parameter values in a
                % formatted way so they can be directly copied into the main script; the
                % -2 means that only fitted parameters are displayed
                disp(' ')
                
                glo2.names = names; % return global to original value
                glo.names_sep = {}; % and put this global back to empty
            end
            
        end
        
        % =================================================================
        % Draw some conclusions from all this fitting
        % =================================================================
        if fit_tox(3) == 3 % only if both controls are requested
            crit = crit_table(nr_fit); % extract criterion that we need
            disp(' ')
            disp('====================================================================')
            for i_set = 1:nd % run through data sets
                
                delL = 2*(mll_coll_cs(i_set) - sum(mll_coll(:,i_set))); % likelihood-ratio criterion that follows chi2
                
                disp(['Conclusions for regular-solvent control, data set: ',num2str(i_set)])
                disp(['Likelihood-ratio criterion (critical value): ',num2str(delL),' (',num2str(crit),')'])
                if delL < crit
                    disp('Conclusion: regular and solvent controls can be lumped')
                    disp('  Note: it is still a good idea to fit the controls with fit_tox(1) = 0.')
                    disp('  This may give a (slightly) different result since the residual standard')
                    disp('  error will be forced equal for both controls.')
                else
                    disp('Conclusion: regular and solvent controls should NOT be lumped')
                end
                if i_set < nd
                    disp('-----------------------------------------------------------')
                end
            end
            disp('====================================================================')
        end
        
        if nd > 1
            disp(' ')
            disp('====================================================================')
            
            for i_c = 1:length(id_ctrl) % run through regular/solvent control
                crit = chi2inv(0.95,nd*nr_fit-nr_fit); % extract criterion that we need
                delL = 2*(mll_coll_d(i_c,1) - sum(mll_coll(i_c,:))); % likelihood-ratio criterion that follows chi2
                if id_ctrl(i_c) == 0
                    str_ctrl = 'regular controls';
                else
                    str_ctrl = 'solvent controls';
                end
                
                disp(['Conclusions for ',str_ctrl,' across all data sets'])
                disp(['Likelihood-ratio criterion (critical value): ',num2str(delL),' (',num2str(crit),')'])
                if delL < crit
                    disp(['Conclusion: ',str_ctrl,' can be lumped across data sets'])
                else
                    disp(['Conclusion: ',str_ctrl,' should NOT be lumped across data sets'])
                end
                if i_c < length(id_ctrl)
                    disp('-----------------------------------------------------------')
                end
            end
            disp('====================================================================')
            disp(' ')
            
            disp('====================================================================')
            for i_c = 1:length(id_ctrl) % run through regular/solvent control
                crit = chi2inv(0.95,nd*nr_fit-nr_fit2); % extract criterion that we need
                delL = 2*(mll_coll_d(i_c,2) - sum(mll_coll(i_c,:))); % likelihood-ratio criterion that follows chi2
                
                str_sep = names_sep{1};
                for i_sep = 2:length(names_sep) % run through extra parameter names
                    str_sep = [str_sep,', ',names_sep{i_sep}];
                end
                
                if id_ctrl(i_c) == 0
                    str_ctrl = 'regular controls';
                else
                    str_ctrl = 'solvent controls';
                end
                disp(['Conclusions for ',str_ctrl,' across all data sets, with individual ',str_sep])
                disp(['Likelihood-ratio criterion (critical value): ',num2str(delL),' (',num2str(crit),')'])
                if delL < crit
                    disp(['Conclusion: ',str_ctrl,' can be lumped across data sets, with separate ',str_sep])
                else
                    disp(['Conclusion: ',str_ctrl,' should NOT be lumped across data sets, with separate ',str_sep])
                end
                if i_c < length(id_ctrl)
                    disp('-----------------------------------------------------------')
                end
            end
            disp('====================================================================')
            
            if fit_tox(3) == 3 % if there are multiple controls ...
                disp(' ')
                disp('====================================================================')
                crit = chi2inv(0.95,nd*nr_fit-nr_fit); % extract criterion that we need
                delL = 2*(mll_coll_d_all - sum(mll_coll_cs2)); % likelihood-ratio criterion that follows chi2
                
                %             disp(' ')
                disp(['Conclusions for lumped controls across all data sets; all parameters common'])
                disp(['Likelihood-ratio criterion (critical value): ',num2str(delL),' (',num2str(crit),')'])
                if delL < crit
                    disp(['Conclusion: both controls can be lumped across all data sets, with one set of parameters'])
                else
                    disp(['Conclusion: both controls should NOT be lumped across all data sets, with one set of parameters'])
                end
                
                crit = chi2inv(0.95,nd*nr_fit-nr_fit2); % extract criterion that we need
                delL = 2*(mll_coll_d_sep - sum(mll_coll_cs2)); % likelihood-ratio criterion that follows chi2
                
                str_sep = names_sep{1};
                for i_sep = 2:length(names_sep) % run through extra parameter names
                    str_sep = [str_sep,', ',names_sep{i_sep}];
                end
                
                disp('-----------------------------------------------------------')
                disp(['Conclusions for lumped controls across all data sets, with individual ',str_sep])
                disp(['Likelihood-ratio criterion (critical value): ',num2str(delL),' (',num2str(crit),')'])
                if delL < crit
                    disp(['Conclusion: both controls can be lumped across all data sets, with separate ',str_sep])
                else
                    disp(['Conclusion: both controls should NOT be lumped across all data sets, with separate ',str_sep])
                end
                disp('====================================================================')
            end
            
            disp(' ')
            disp('NOTE: these conclusions are based on statistics only (likelihood-ratio tests).')
            disp('A statistically significant difference is not necessarily biologically meaningful.')
            disp('Therefore, also pay close attention to the visual fits.')
            disp('And note that the survival data are included in these fits/comparisons.')
            
        end
        disp(' ')
        
    case -1 % for fitting the survival control response only ...
        % This piece of code selects only the treatment with identifier 0
        % (zero), removes all data apart from survival, and turns fitting
        % off for all parameters but background hazard.
        
        if par.hb(2) == 0
            disp('Background hazard is set to non-fitted, so we can stop immediately!')
            par_out = par;
            return
        end
        
        % Never show iterations for these control fits
        opt_optim.it = 0; % show iterations of the simplex optimisation (1, default) or not (0)
        
        % The exact setting for moa and feedb is not relevant for the
        % control, but they must be defined to prevent errors.
        glo.moa   = MOA(1,:);   % use FIRST entry in MOA matrix
        glo.feedb = FEEDB(1,:); % use FIRST entry in FEEDB matrix
        nd        = size(DATA,1); % number of data sets
        
        loc_all = [];
        for i_set = 1:nd % run through data sets
            [~,loc_i] = ismember(id_ctrl+(i_set-1)*100,X0mat(1,:)); % find where controls are in X0mat in this data set
            loc_all = cat(2,loc_all,loc_i); % add to previous locations
        end
        loc_all(loc_all == 0) = []; % remove zero from location vector (e.g., when there is no solvent control)
        X0mat = X0mat(:,loc_all); % use only those columns for X0mat
        
        for i_sets = 1:size(DATA,1)
            DATA(i_sets,[glo.locD glo.locL glo.locR]) = {[0] [0] [0]}; % remove all but the survival data from the complete data set
            % Note: this also works in case there are multiple data sets per state. 
        end
        
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
            par.hb(1) = 0.01; % and start from a reasonable default (if parspace explorer is NOT used)
            if isfield(par,'a') && par_rem.a(2) == 1 % assume a is for a Weibull background coefficient!
                par.a(2) = 1; % fit that coefficient as well
            end
        end
        
        glo2.ctot = X0mat(1,:); % make sure that the total concentration vector to analyse only has control identifiers
        
        opt_plot.statsup = [glo.locD glo.locL glo.locR]; % vector with states to suppress in plotting fits
        % No need for annotations or legends as we fit hb only
        opt_plot.annot  = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
        opt_plot.legsup = 0; % set to 1 to suppress legends on fits
        glo.basenm      = [basenm_rem,'_hb']; % modify basename to make a separate PS.mat file
        
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
        
    case 0 % for fitting the control response only ...
        % This piece of code selects only the treatment with identifier 0
        % (zero), and turns fitting off for the toxicity parameters.
        
        % Never show iterations for these control fits
        opt_optim.it = 0; % show iterations of the simplex optimisation (1, default) or not (0)
        
        % The exact setting for moa and feedb is not relevant for the
        % control, but they must be defined to prevent errors.
        glo.moa   = MOA(1,:);   % use FIRST entry in MOA matrix
        glo.feedb = FEEDB(1,:); % use FIRST entry in FEEDB matrix
        nd        = size(DATA,1); % number of data sets
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            if opt_optim.type ~= 4 % for simplex fitting
                opt_optim.fit = 0; % fit the parameters (1), or don't (0)
            else
                opt_optim.ps_saved = 1; % use saved set for parameter-space explorer (1) or not (0);
                opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
                % Fit needs to be 1 otherwise the saved set is not loaded.
            end
        end
        
        loc_all = [];
        for i_set = 1:nd % run through data sets
            [~,loc_i] = ismember(id_ctrl+(i_set-1)*100,X0mat(1,:)); % find where controls are in X0mat in this data set
            loc_all = cat(2,loc_all,loc_i); % add to previous locations
        end
        loc_all(loc_all == 0) = []; % remove zero from location vector (e.g., when there is no solvent control)
        X0mat = X0mat(:,loc_all); % use only those columns for X0mat
        
        pmat      = packunpack(1,par,[],WRAP); % turn parameter structure into a matrix
        pmat(ind_tox:end,2) = 0; % don't fit the toxicity parameters now
        par       = packunpack(2,[],pmat,WRAP); % turn parameter matrix back into a structure
        glo2.ctot = X0mat(1,:); % make sure that the total concentration vector to analyse only has control identifiers
        
        opt_plot.statsup = [glo.locD glo.locS]; % vector with states to suppress in plotting fits
        % No need for annotations or legends as we fit basic parameters only
        opt_plot.annot  = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
        opt_plot.legsup = 0; % set to 1 to suppress legends on fits
        
        glo.basenm      = [basenm_rem,'_ctrl']; % modify basename to make a separate PS.mat file
        
        % assume that hb has been fitted separately ... don't do it again
        par.hb(2)        = 0; % don't fit the background hazard rate
        if isfield(par,'a') % assume a is for a Weibull background coefficient!
            par.a(2) = 0; % don't fit that coefficient as well
        end
        for i_sets = 1:size(DATA,1)
            DATA(i_sets,glo.locS) = {[0]}; % remove the survival data from the complete data set
            % this should also work for multiple sets per state
        end
        
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
        
    case 1 % for fitting tox parameters only ...
        
        % opt_plot.statsup = [glo.locD]; % vector with states to suppress in plotting fits
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            if opt_optim.type ~= 4 % for simplex fitting
                opt_optim.fit = 0; % fit the parameters (1), or don't (0)
            else
                opt_optim.ps_saved = 1; % use saved set for parameter-space explorer (1) or not (0);
                opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
                % Fit needs to be 1 otherwise the saved set is not loaded.
            end
        else
            % This piece of code turns fitting off for the basic parameters so only
            % the toxicity parameters are fitted.
            pmat = packunpack(1,par,[],WRAP); % turn output parameter structure into a matrix
            pmat(1:ind_tox-1,2) = 0; % turn fitting for the basic parameters off in this vector
            par = packunpack(2,[],pmat,WRAP); % turn parameter matrix into a structure
        end
        
        nd = size(DATA,1); % number of data sets
        loc_ctrl_all = [];
        loc_solv_all = [];
        for i_set = 1:nd % run through data sets
            [~,loc_ctrl] = ismember(0+(i_set-1)*100,X0mat(1,:)); % find where controls are in X0mat in this data set
            [~,loc_solv] = ismember(id_solvent+(i_set-1)*100,X0mat(1,:)); % find where solvent controls are in X0mat in this data set
            loc_ctrl_all = cat(2,loc_ctrl_all,loc_ctrl); % add to previous locations
            loc_solv_all = cat(2,loc_solv_all,loc_solv); % add to previous locations
        end
        % remove indices that are zero!
        loc_ctrl_all(loc_solv_all==0) = [];
        loc_solv_all(loc_solv_all==0) = [];
        
        switch fit_tox(3) 
            case 1 % only use regular control, so remove solvent control
                % X0mat(:,X0mat(1,:)==loc_solv_all) = []; % remove solvent control(s)
                if sum(loc_solv_all) > 0
                    X0mat(:,loc_solv_all) = []; % remove solvent control(s)
                end
                glo2.ctot = X0mat(1,:); % update the concentration vector that will be used!
            case 2 % only use solvent control, so remove regular control
                % X0mat(:,X0mat(1,:)==loc_ctrl_all) = []; % remove regular control(s)
                if sum(loc_ctrl_all) > 0
                    X0mat(:,loc_ctrl_all) = []; % remove regular control(s)
                end
                glo2.ctot = X0mat(1,:); % update the concentration vector that will be used!
        end
        % Note: controls are never used for fitting at this stage, but now
        % the unused control is also removed from the plots.
        
        % Run through the permutations of MoA and feedbacks. The basename,
        % as used for naming the mat file and the output plots, is adapted
        % to identify the settings.
        MLL_coll = zeros(size(MOA,1),size(FEEDB,1)); % initialise matrix to catch MLLs
        par_coll = cell(size(MOA,1),size(FEEDB,1));  % initialise matrix to catch parameter sets
        for i = 1:size(MOA,1) % run through all MoAs
            for j = 1:size(FEEDB,1) % run through all feedback configuration
                glo.moa    = MOA(i,:); % change global for MoA
                glo.feedb  = FEEDB(j,:); % change global for feedback configuration
                glo.basenm = [basenm_rem,'_moa',sprintf('%d',glo.moa),'_feedb',sprintf('%d',glo.feedb)];
                
                % automatic creation of search ranges for DEBtox tox parameters
                if fit_tox(2) == 1 && opt_optim.type == 4 && skip_sg == 0 % only when using parspace explorer for fitting
                    par = startgrid_debtox(par); % create search ranges for the tox parameters
                end % this needs to be in the loop as the MoA affects the ranges!
                
                % optimise and plot (fitted parameters in par_out)
                [par_tmp,MLL] = calc_optim(par,opt_optim); % start the optimisation
                
                names      = fieldnames(par); % extract all field names of par
                names_tmp  = fieldnames(par_tmp); % extract all field names of par_out from MAT file
                ind_fittag = ~strcmp(names_tmp,'tag_fitted');
                names_tmp  = names_tmp(ind_fittag); % make sure that the fit tag is not in names_out
                ind_fittag = ~strcmp(names,'tag_fitted');
                names      = names(ind_fittag); % make sure that the fit tag is not in names
                
                if isequal(names,names_tmp)
                    calc_and_plot(par_tmp,opt_plot); % calculate model lines and plot them
                else
                    warning('No plot is produced in automatic_runs since the parameters in the saved MAT file do not match those in the workspace')
                    % This may happen for validation, where a saved MAT
                    % file is not exactly the same as the current set.
                    % Plotting par_out with calc_and_plot would then lead
                    % to unexpected output, or even errors.
                end
                MLL_coll(i,j) = MLL; % collect the MLL so we can compare them
                par_coll{i,j} = par_tmp; % collect the par so we can output the best one
                drawnow
            end
        end
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            AIC_coll = nan(size(MOA,1),size(FEEDB,1)); % use NaNs matrix to catch MLLs
            % If you did not fit, AIC is meaningless
            if size(MOA,1) == 1 && size(FEEDB,1)
                best_MoaFb = [1 1]; % indices for best configuration
                par_out    = par_coll{1,1}; % output the best par
            end
            % If a script is run to display a saved set, we can safely say
            % it's the best (that may prevent an unnecessary error).
        else
            AIC_coll  = 2*sum(pmat(:,2))+2*MLL_coll; % calculate AIC from MLL
            AIC_coll  = AIC_coll - min(AIC_coll(:)); % calculate delta AIC relative to best one
        end
        loss_coll = exp(-AIC_coll/2); % relative probability to minimise information loss
        
        % display on screen and in diary
        diary ('results.out') % collect output in the diary "results.out"
        disp(' ')
        disp('MoA  feedbacks   MLL     delta-AIC     prob.    best')
        disp('====================================================================')
        for i = 1:size(MOA,1) % run through all MoAs
            for j = 1:size(FEEDB,1) % run through all feedback configuration
                str1 = sprintf('%d',MOA(i,:)); % string for the MoA configuration
                str2 = sprintf('%d',FEEDB(j,:)); % string for the feedbacks configuration
                fprintf('%s  %s %#10.2f %#10.2f %#10.2f ',str1,str2,MLL_coll(i,j),AIC_coll(i,j),loss_coll(i,j))
                if AIC_coll(i,j) == 0 % if this is the best one ...
                    fprintf('%s','     *') % give it a star
                    best_MoaFb = [i j]; % indices for best configuration
                    par_out    = par_coll{i,j}; % output the best par
                end
                fprintf('\n')
            end
            disp('====================================================================')
        end
        diary off  % turn off the diary function

end

% Return settings to their original value
if fit_tox(1) ~= 1 % when we fit tox data, it is good to keep the modified glo
    glo   = glo_rem;
end
glo2      = glo2_rem;
X0mat     = X0mat_rem;
DATA      = DATA_rem;
W         = W_rem;
drawnow