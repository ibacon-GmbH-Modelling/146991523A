function par_best = start_vals_guts(par)

% This function searches for starting values for the simple GUTS cases in
% the GUTS package (only reduced model for constant exposure). Using "rules
% of thumb", a range of values is tried (with a rough optimisation), and
% the best one displayed. It only uses the survival data set indicated by
% glo.locS (and if there are multiple, only the first one).
% 
% Note that the openGUTS software (Matlab-based and standalone) is better
% suited to find the global optimum (as well as the parameter-space
% explorer as implemented in BYOM). Nevertheless, this function will
% generally be good enough for well-behaved cases, and could be extended in
% the future.
% 
% Author : Tjalling Jager
% Date: May 2020
% Web support: http://www.debtox.info/byom.html
 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global DATA glo

opttype = 2; % use experimental swarm optimisation (1) or multiple Simplex runs (2)

pmat   = packunpack(1,par,0);  % transform structure into a regular matrix

% From the data set, extract the time and concentration vectors from the
% first data set (survival) only.
C_data = DATA{1,glo.locS}(1,2:end); % concentration vector in the survival data set
T_data = DATA{1,glo.locS}(2:end,1); % time vector in the survival data set
S_data = DATA{1,glo.locS}(2:end,2:end); % survivors

% For each parameter that is fitted, come up with a reasonable range to try
if par.mw(2) == 1
    mw_try = mean([C_data(1:end-1);C_data(2:end)]); % for threshold, try values between the exposure concentrations
    mw_try = mw_try(S_data(end,1:end-1)./S_data(1,1:end-1)>0.5); % only keep the values were there is less than 50% effect at the end
    mw_try = max(par.mw(3),mw_try);
    mw_try = min(par.mw(4),mw_try);
    mw_try = unique(mw_try); % make sure that min-max does not lead to same values
    switch length(mw_try) % make sure that mw_try is at least 3 values long
        case 0, mw_try = linspace(par.mw(3),min(C_data(C_data>0)),3);
        case 1, mw_try = linspace(par.mw(3),mw_try,3);
        case 2, mw_try = [mw_try(1) mean(mw_try) mw_try(2)];
    end
else
    mw_try = par.mw(1);
end
if par.kd(2) == 1
    Fss = [0.95 0.7 0.5 0.2 0.05]; % fraction of steady state within test duration
    kd_try = -max(T_data)*log(1-Fss); % translate to ke values
    kd_try = max(par.kd(3),kd_try); % make sure the test range is within min/max bounds
    kd_try = min(par.kd(4),kd_try);
    kd_try = unique(kd_try); % make sure that min-max does not lead to same values
else % if it is not fitted, only 'try' the given value in the byom script
    kd_try = par.kd(1);
end
if par.bw(2) == 1
    bw_try = [1 0.33 0.1 0.033 0.01] * 10/max(T_data) * 10/max(C_data); % rule of thumb for reasonable killing rates
    bw_try = max(par.bw(3),bw_try);
    bw_try = min(par.bw(4),bw_try);
    bw_try = unique(bw_try); % make sure that min-max does not lead to same values
else
    bw_try = par.bw(1);
end
if par.Fs(2) == 1
    Fs_try = [1.2 2 4 6]; % reasonable range for fraction spread
    Fs_try = max(par.Fs(3),Fs_try);
    Fs_try = min(par.Fs(4),Fs_try);
    Fs_try = unique(Fs_try); % make sure that min-max does not lead to same values
else
    Fs_try = par.Fs(1);
end
if par.hb(2) == 1
    hb_try = [0 0.01]; % reasonable range for background mortality
    hb_try = max(par.hb(3),hb_try);
    hb_try = min(par.hb(4),hb_try);
    hb_try = unique(hb_try); % make sure that min-max does not lead to same values
else
    hb_try = par.hb(1);
end

switch opttype
    
    case 1
        
        %% Swarm optimisation
        
        mw_range = [min(mw_try) max(mw_try)];
        kd_range = [min(kd_try) max(kd_try)];
        bw_range = [min(bw_try) max(bw_try)];
        Fs_range = [1.2 6];
        hb_range = [0.001 0.02];
        
        par_tmp = par; % remember the parameter structure, so we can mess with a copy
        % put the estimated ranges in, and put all on normal scale (it does not
        % seem to like negative values ...
        par_tmp.mw = [mean(mw_range) par.mw(2) mw_range 1];
        par_tmp.kd = [mean(kd_range) par.kd(2) kd_range 1];
        par_tmp.bw = [mean(bw_range) par.bw(2) bw_range 1];
        par_tmp.Fs = [mean(Fs_range) par.Fs(2) Fs_range 1];
        par_tmp.hb = [mean(hb_range) par.hb(2) hb_range 1];
        
        pmat_tmp   = packunpack(1,par_tmp,0);  % transform structure into a regular matrix
        
        % don't put parameters on log scale
        pfit   = pmat_tmp(pmat_tmp(:,2)==1,1); % parameter values that are to be estimated
        % set parameters based on dimensionality of the problem, and width of min-max bounds!
        nobirds = 200; % number of birds in the swarm
        noiter  = 20; % number of iterations (don't need too many as Simplex will follow)
        swarm_bnds = pmat_tmp(:,[3 4]);
        swarm_bnds(pmat_tmp(:,5)==0,:) = log10(swarm_bnds(pmat_tmp(:,5)==0,:)); % if parameters are estimated on log-scale, the bounds must also be log scale
        swarm_bnds = swarm_bnds(pmat_tmp(:,2)==1,:); % only fitted parameters
        phat = Particle_Swarm_Optimization (nobirds,length(pfit),swarm_bnds,@transfer,'min',2,2,2,0.4,0.9,noiter,pmat_tmp);
        
        pmat_tmp(pmat_tmp(:,2)==1,1) = phat; % put best estimate back (all on normal scale
        pmat_tmp = [pmat_tmp(:,1) pmat(:,[2:5])];
        par_best = packunpack(2,0,pmat_tmp); % transform matrix into a structure to return to the script
        
    case 2
        
        %% Multiple Simplex optimisations
        
        oldopts    = optimset('fminsearch'); % default options for fminsearch
        options_s1 = optimset(oldopts,'FunValCheck','on'); % function check on
        options_s1 = optimset(options_s1,'Display','off','FunValCheck','on'); % do not show iterations
        options_s2 = optimset(options_s1,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50); % rough Simplex options
        
        % disp(' ')
        % disp('Finding starting values for the parameters ...')
        
        f = waitbar(0,'Finding starting values for the GUTS parameters. Please wait.','Name','start_vals_guts.m');
        trials_tot = length(mw_try)*length(kd_try)*length(bw_try)*length(Fs_try)*length(hb_try);        
        
        teller = 1;
        for i = 1:length(mw_try)
            for j = 1:length(kd_try)
                for k = 1:length(bw_try)
                    for m = 1:length(Fs_try)
                        for n = 1:length(hb_try)
                            par_tmp = par; % remember the parameter structure, so we can mess with a copy
                            par_tmp.mw(1) = mw_try(i);
                            par_tmp.kd(1) = kd_try(j);
                            par_tmp.bw(1) = bw_try(k);
                            par_tmp.Fs(1) = Fs_try(m);
                            par_tmp.hb(1) = hb_try(n);
                            
                            pmat   = packunpack(1,par_tmp,0);  % transform structure into a regular matrix
                            % put parameters that need to be on log scale on log10 scale
                            pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));
                            pfit   = pmat(pmat(:,2)==1,1); % parameter values that are to be estimated
                            
                            [phat,FVAL,EXITFLAG,OUTPUT] = fminsearch('transfer',pfit,options_s2,pmat); % rough estimation
                            COLL(teller,1) = FVAL; % remember the minloglik value
                            COLL(teller,2:length(phat)+1) = phat'; % and remember the parameter vector
                            teller = teller+1;
                            
                            waitbar(teller/trials_tot,f) % make a nice waiting bar
                            
                        end
                    end
                end
            end
        end
        
        close(f) % close the waiting bar
        
        ind_best = find(COLL(:,1)==min(COLL(:,1)),1,'first'); % where is the best value?
        phat = (COLL(ind_best,2:end))'; % extract that parameter vector
        
        % minimisation is done with the reduced parameter set, other parameters in par are kept
        % to the fixed value. This means that phat will also contain the limited set only.
        % To make the vector whole again:
        p = pmat(:,1);          % take the first column with the initial parameter values
        p(pmat(:,2)==1) = phat ; % replace values in total vector with values that were fitted
        
        % put parameters that need to be fitted on log scale back on normal scale
        p(pmat(:,5)==0)  = 10.^(p(pmat(:,5)==0));
        
        p = max(p,pmat(:,3));    % force values within the specified lower parameter bounds
        p = min(p,pmat(:,4));    % force values within the specified upper parameter bounds
        pmat(:,1) = p;           % return the new p vector into the matrix pmat
        
        par_best = packunpack(2,0,pmat); % transform matrix into a structure to return to the script
        
end

% disp(' ')
% disp(['Best parameter set from the starting value estimator, minloglik: ',num2str(COLL(ind_best,1))])
% format short g
% disp(par_best)

