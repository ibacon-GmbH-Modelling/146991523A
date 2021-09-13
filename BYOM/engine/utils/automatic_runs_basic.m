function par_out = automatic_runs_basic(fit_tox,par,ind_tox,skip_sg,opt_optim,opt_plot,varargin)

% This is a helpful function to run a TKTD model analysis in parts, based
% on the setting in fit_tox. We can fit the basic (non-tox) parameters on
% the control, or the tox parameters only to the complete data set. Note
% that this code expects the control to have identifier zero.
% 
% In general, this function will be used in conjunction with the
% parameter-space explorer (opt_optim.type = 4), at least for fitting the
% toxicity data, but most of it will also work with other optimisation
% routines.
% 
% At this moment, this code is for GUTS, GUTS-immobility and DEBtox2019
% analyses only. I may adapt it in the future for other analyses as well.
% Note that GUTS analysis is triggered by entering ind_tox=-1, and
% GUTS-immobility by ind_tox=-2. This is possible since there is only 1
% basic parameter (hb).
% 
% Select whether to fit the control treatment (0), the toxicity treatments
% (1; the control is shown as well), or both together (2). A good strategy
% is to fit the basic parameters to the control data first (fit_tox=0),
% copy the best values into the parameter matrix above, and then fit the
% toxicity parameter to the complete data set (fit_tox=1). The code below
% automatically keeps the parameters fixed that need to be fixed.
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global X0mat glo glo2 DATA

% To ensure backwards compatability (since I added the MoA 'hazards to repro' again)
if isfield(glo,'moa') && length(glo.moa) == 4
    glo.moa(5) = 0; % make sure the length of MOA is 5
end

% remember the globals as we will change them
X0mat_rem = X0mat;
glo2_rem  = glo2;
DATA_rem  = DATA;

switch fit_tox
    case 0 % for fitting the control response only ...
        % This piece of code selects only the treatment with identifier 0
        % (zero), and turns fitting off for the toxicity parameters.
        X0mat = X0mat(:,X0mat(1,:)==0); % only keep controls
        pmat  = packunpack(1,par,[]); % turn parameter structure into a matrix
        
        switch ind_tox
            case -1 % this signals that we have a GUTS analysis
                pmat(:,2) = 0; % don't fit ANY parameters now
                par       = packunpack(2,[],pmat); % turn parameter matrix back into a structure
                par.hb(2) = 1; % but DO fit the background hazard rate
            case -2 % this signals that we have a GUTS-immobility analysis
                % Specifically for the immobility package: remove all data other
                % than those for the dead animals. However, startgrid_guts will
                % only use the time vector of the first one for its derivation of
                % search ranges.
                for i_sets = 1:size(DATA,1)
                    DATA(i_sets,[glo.locC glo.locD glo.loc_h glo.loc_i glo.loc_id]) = {[0] [0] [0] [0] [0]}; % remove all but the survival data from the complete data set
                    % Note: this also works in case there are multiple data sets per state.
                end
                opt_plot.statsup = [glo.locC glo.locD glo.loc_h glo.loc_i glo.loc_id]; % vector with states to suppress in plotting fits
                pmat(:,2) = 0; % don't fit ANY parameters now
                par       = packunpack(2,[],pmat); % turn parameter matrix back into a structure
                par.hb(2) = 1; % but DO fit the background hazard rate
            otherwise % we have a DEBtox analysis
                pmat(ind_tox:end,2) = 0; % don't fit the toxicity parameters now
                par   = packunpack(2,[],pmat); % turn parameter matrix back into a structure
        end
        
        glo2.ctot = 0; % make sure that the total concentration vector to analyse only has a zero
        % No need to update WRAP, since packunpack only needs names
        
        % optimise and plot (fitted parameters in par_out)
        par_out = calc_optim(par,opt_optim); % start the optimisation
        calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

        print_par(par_out,-2) % this prints out the optimised parameter values in a
        % formatted way so they can be directly copied into the code above; the
        % -2 means that only fitted parameters are displayed
        
        if opt_optim.type ~= 4 && ~isempty(varargin)
            opt_prof = varargin{1};
            calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
        end
        
    case 1 % for fitting tox parameters only ...
        % This piece of code turns fitting off for the basic parameters so only
        % the toxicity parameters are fitted.
        
        switch ind_tox
            case -1 % this signals that we have a GUTS analysis
                par.hb(2) = 0; % turn fitting for the background hazard off
                if opt_optim.type == 4 && skip_sg ~= 1 % for the parspace explorer ...
                    par = startgrid_guts(par); % replace par by estimates based on data set (as used by openGUTS)
                end
            case -2 % this signals that we have a GUTS-immobility analysis
                par.hb(2) = 0; % turn fitting for the background hazard off
                if opt_optim.type == 4 && skip_sg ~= 1 % for the parspace explorer ...
                    par = startgrid_guts_immob(par); % replace par by estimates based on data set (as used by openGUTS)
                end
            otherwise % we have a DEBtox analysis
                pmat = packunpack(1,par,[]); % turn output parameter structure into a matrix
                pmat(1:ind_tox-1,2) = 0; % turn fitting for the basic parameters off in this vector
                par  = packunpack(2,[],pmat); % turn parameter matrix into a structure
                if opt_optim.type == 4 && skip_sg ~= 1 % for the parspace explorer ...
                    par = startgrid_debtox(par); % create search ranges for the fitted tox parameters
                end
        end
                
        % optimise and plot (fitted parameters in par_out)
        par_out = calc_optim(par,opt_optim); % start the optimisation
        calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
        
    case 2 % for fitting all parameters
        if opt_optim.type == 4 && skip_sg ~= 1 % for the parspace explorer ...
            switch ind_tox
                case -1 % this signals that we have a GUTS analysis
                    par = startgrid_guts(par); % replace par by estimates based on data set (as used by openGUTS)
                case -2 % this signals that we have a GUTS-immobility analysis
                    par = startgrid_guts_immob(par); % replace par by estimates based on data set (as used by openGUTS)
                otherwise % we have a DEBtox analysis
                    par = startgrid_debtox(par); % create search ranges for the fitted tox parameters
            end
        end  
        
        % optimise and plot (fitted parameters in par_out)
        par_out = calc_optim(par,opt_optim); % start the optimisation
        calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
end

% return the globals to their original values
X0mat = X0mat_rem;
glo2  = glo2_rem;
DATA  = DATA_rem;
drawnow
