function [par_out,FVAL,GoodFit,pmat_print] = calc_optim(par,opt_optim)

% Usage: [par_out,FVAL,GoodFit,pmat_print] = calc_optim(par,opt_optim)
%
% This function starts up the optimisation of the parameter values; it is
% called from <calc_and_plot.m>. It calls the optimisation function as
% specified in the options structure. By default: the Matlab function
% <fminsearch>, which performs a Nelder-Mead Simplex search. This file also
% display the results in a formatted manner. Optional outputs are the
% best-fitting parameters, the min-log-likelihood and the r-square.
%
% The structure <opt_optim> can be used to pass options to the optimisation
% routine (see <prelim_checks.m>).
%
% The parameter structure <par> is input, output is the structure <par_out>
% which contains the optimised values in the first position, and <FVAL>
% which contains the minus log-likelihood value.
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global DATA glo glo2

% extract options from the opt_optim provided (filled in prelim_checks)
fit       = opt_optim.fit;     % fit the parameters (1), or don't (0)
it        = opt_optim.it;      % show iterations of the optimisation (1, default) or not (0)
opttype   = opt_optim.type;    % optimisation method 1) simplex, 2) simulated annealing, 3) swarm optimisation
simno     = opt_optim.simno;   % for simplex: number of runs of fminsearch (starting again from previous best)
anntemp   = opt_optim.anntemp; % for simulated annealing, stopping temperature
swno      = opt_optim.swno;    % for swarm optimisation, number of particles
swit      = opt_optim.swit;    % for swarm optimisation, number of iterations

filenm     = glo.basenm;        % name for indentification of the analysis
pmat_print = []; % this output is normally returned as empty, but when parspace is used, it contains the CIs

% extract parameters from the general globals, glo and glo2
n_D       = glo2.n_D;
n_X       = glo2.n_X;
names     = glo2.names;
namesz    = glo2.namesz;

pmat = packunpack(1,par,0);  % transform structure into a regular matrix
ind_log = find(pmat(:,5)==0); % find parameters that are on log scale
pmat(ind_log,1) = log10(pmat(ind_log,1)); % put parameters that need to be on log scale on log10 scale
pfit = pmat(pmat(:,2)==1,1); % parameter values that are to be estimated (log10 scale where needed)

if isempty(pfit) || fit == 0 % if there is nothing to be fitted anyway ...
    if opttype == 4 && opt_optim.ps_saved == 1
        % then there is no problem to NOT fit any parameters: this function
        % will be used to load a previously saved sample from parameter
        % space (and to display the results as well).
    else
        par_out    = par; % simply return same matrix (but without fitted tag)
        % In principle, it would be possible to use load_rnd to return par
        % from a saved mat file, but that requires opt_conf.
        FVAL       = -1;
        GoodFit    = [];
        pmat_print = [];
        disp(['Time required: ' secs2hms(toc)])
        return % go back to the main script ...
    end
end

oldopts    = optimset('fminsearch'); % default options for fminsearch
options_s1 = optimset(oldopts,'FunValCheck','on'); % function check on
if it == 1 % show iterations
    options_s1 = optimset(options_s1,'Display','iter','FunValCheck','on'); % show all iterations
else
    options_s1 = optimset(options_s1,'Display','off','FunValCheck','on'); % show no output at all
end
options_s2 = optimset(options_s1,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',30*length(pfit)); % rough Simplex options

% Call the optimisation routine as specified by opttype
switch opttype
    case 1 % standard, two (or more) consecutive simplexes
        
        [phat,FVAL,EXITFLAG,OUTPUT] = fminsearch('transfer',pfit,options_s2,pmat); % first rough estimation
        out_iter = OUTPUT.iterations; % remember iterations
        % disp(['rough fminsearch iteration:',num2str(OUTPUT.iterations)])
        
        for no_opt = 1:simno-1        
            [phat,FVAL,EXITFLAG,OUTPUT] = fminsearch('transfer',phat,options_s1,pmat); % then a more detailed one
            out_iter = out_iter+OUTPUT.iterations; % count total iterations
        end
        % disp(['std fminsearch iteration:',num2str(OUTPUT.iterations)])
        
    case 2 % simulated annealing + simplex (experimental)
        
        opt_ann.Verbosity = 2; % controls output on screen
        opt_ann.StopTemp  = anntemp; % crude temperature for stopping as Simplex will finish it off
        [phat,~] = anneal(@transfer,pfit,opt_ann,pmat); % first rough estimation
        [phat,FVAL,EXITFLAG,OUTPUT]  = fminsearch('transfer',phat,options_s1,pmat); % followed by a detailed simplex
        out_iter = OUTPUT.iterations; % remember iterations
        
    case 3 % particle swarm + simplex (experimental)
        
        % set parameters based on dimensionality of the problem, and width of min-max bounds!
        nobirds = swno; % number of birds in the swarm
        noiter  = swit; % number of iterations (don't need too many as Simplex will follow)
        swarm_bnds = pmat(:,[3 4]); 
        swarm_bnds(ind_log,:) = log10(swarm_bnds(ind_log,:)); % if parameters are estimated on log-scale, the bounds must also be log scale
        swarm_bnds = swarm_bnds(pmat(:,2)==1,:); % only fitted parameters
        phat = Particle_Swarm_Optimization (nobirds,length(pfit),swarm_bnds,@transfer,'min',2,2,2,0.4,0.9,noiter,pmat);
        out_iter = 0; % do not count swarm iterations in noiter
        
        [phat,FVAL,EXITFLAG,OUTPUT]  = fminsearch('transfer',phat,options_s1,pmat); % followed by a detailed simplex
        out_iter = out_iter+OUTPUT.iterations; % remember iterations
         
    case 4 % Parameter-space explorer
        % This is the parameter-space explorer (as also used for the
        % openGUTS software) in a more general setting. It uses Matlab
        % functions in a separate sub-folder <parspace> of the engine.
        
        [pmat,pmat_print,FVAL,names_out] = calc_optim_ps(pmat,opt_optim);
        % NOTE PMAT IS ON LOG-SCALE WHERE NEEDED AND PMAT_PRINT ON NORMAL SCALE
        phat = pmat(pmat(:,2)==1,1); % extract fitted parameters from pmat
        
        EXITFLAG = 1; % this can be sorted out later, but for now assume it always exits okay
        out_iter = inf; % iterations are not really interesting for PS explorer
        
        ind_log = find(pmat(:,5)==0); % find parameters that are on log scale
        % this needs to be done again, as calc_optim_ps may modify the
        % log-settings of the threshold when slow kinetics is indicated!
               
        if opt_optim.ps_saved == 1 % when using a saved set ...
            % create pfit from saved set rather than input set (this allows recontrsuction of the original AIC)
            pfit = pmat(pmat(:,2)==1,1); % parameter values that are to be estimated (log10 scale where needed)
        end
        
    case 5 % Alternative Nelder-Mead simplex
        % This implements an alternative to Matlab's fminsearch, namely the
        % function nelmin (see code of this function for more information).
        % The reason to include this option is that this alternative is
        % used in the openGUTS software, so it can be tested in BYOM as
        % well.
        
        [phat,FVAL,icount,numres,ifault] = nelmin(@transfer,length(pfit),pfit,1e-8,0.05*pfit,1,30*length(pfit),pmat);
        out_iter = icount; % remember iterations
        % disp(['rough nelmin iteration:',num2str(icount)])
        
        for no_opt = 1:simno-1        
            [phat,FVAL,icount,numres,ifault] = nelmin(@transfer,length(pfit),pfit,1e-15,0.05*pfit,1,200*length(pfit),pmat);
        end
        out_iter = out_iter+icount; % count total iterations
        % disp(['rough nelmin iteration:',num2str(icount)])
        
        % make sure the exitflag is defined in the same manner as with fminsearch
        if ifault == 0
            EXITFLAG = 1;
        else
            EXITFLAG = 0;
        end
end

if opttype == 4 && opt_optim.ps_saved == 1 && ~isempty(names_out)
    % if this series of arguments is true, we are likely called to display
    % results from a saved set, in a validation setting (with a slightly
    % different model than the one used for calibration). calc_optim_ps
    % returns the parameters from the saved set ...
    GoodFit = NaN;
    zvd     = [];
    names   = names_out;
    glo2_rem   = glo2;
    glo2.names = names; % also change the glo2 ...
else
    glo2_rem = [];
    [~,GoodFit,zvd] = transfer(phat,pmat); % get goodness of fit measures (and zvd output, which is used in an eval below)
    disp('')
    disp('Goodness-of-fit measures for each state and each data set (R-square)')
    disp('Based on the individual replicates, accounting for transformations, and including t=0')
    % disp(GoodFit)
    find_data = reshape(1:length(GoodFit),n_D,n_X); % helper to find which state a data set is
    for i = 1:n_X % run through states
        fprintf('state %#2.0f',i)
        for j = 1:n_D % run through number of data sets per state
            fprintf(' %#6.2f',GoodFit(find_data(j,i)))
        end
        fprintf('\n')
    end
end

survdata = 0; % check if there are survival data
cntdata  = 0; % counter for real data sets (bigger than 1 element)
for i = 1:n_X % run through states
    for j = 1:n_D % run through number of data sets per state
        if DATA{j,i}(1,1) < 0 % likelihood type
            survdata = 1;
        end
        if prod(size(DATA{j,i})-1) > 1
            cntdata = cntdata + 1;
        end
    end
end 

if survdata == 1
    warning('off','backtrace') % no need to display where the warning is generated
    warning('R-square is not appropriate for survival data so interpret qualitatively!')
    warning('(R-square not calculated when missing/removed animals; use plot_tktd)')
    warning('on','backtrace'), disp(' ')
end
% GoodFit is a vector with goodness-of-fit measures for each data set. For
% all data sets, it is now the R-squared (non-adjusted). For survival data
% it is not very appropriate, but may serve to give a rough idea of how
% well the survival probability corresponds to the observed frequencies.

% minimisation is done with the reduced parameter set, other parameters in
% par are kept to the fixed value. This means that phat will also contain
% the limited set only. To make the vector whole again:
p0 = pmat(:,1);          % take the first column with the *initial* parameter values
p  = p0;                 % copy the initial values to a new vector
p(pmat(:,2)==1) = phat ; % replace values in total vector with values that were fitted

if opttype == 4 || opttype == 3
    % for these optimisation routines, the initial values are meaningless
    p0(:) = nan; % so just represent them as NaN
end    

% put parameters that need to be fitted on log scale back on normal scale
p(ind_log)  = 10.^p(ind_log);
p0(ind_log) = 10.^p0(ind_log);

p = max(p,pmat(:,3));    % force values within the specified lower parameter bounds
p = min(p,pmat(:,4));    % force values within the specified upper parameter bounds
pmat(:,1) = p;           % return the new p vector into the matrix pmat

par_out = packunpack(2,0,pmat); % transform matrix into a structure to return to the script
par_out.tag_fitted = 'y';       % add a tag to let plotting function know it is fitted
if ~isempty(glo2_rem)
    glo2 = glo2_rem;
end

% ========== Now display the results on screen ==========
diary (glo.diary) % collect output in the diary "results.out"
% disp(' ')
disp('=================================================================================')
disp('Results of the parameter estimation with BYOM version 6.0 beta 4')
disp('=================================================================================')
fprintf('   Filename      : %s \n',filenm) % name of the mydata file
ti = clock; % what is the current time?
fprintf('   Analysis date : %s (%02.0f:%02.0f) \n',date,ti(4),ti(5));

disp('   Data entered  :')
for i = 1:n_X % run through states
    for j = 1:n_D % run through number of data sets per state
        lam    = DATA{j,i}(1,1); % likelihood type
        datsiz = size(DATA{j,i})-1; % number of rows-columns of observations
        if n_D == 1 % if there is only one data set ...
            fprintf('     data state %1.0f: %1.0fx%1.0f, ',i,datsiz(1),datsiz(2))
        else % if there's more, also print the number of the data set between parentheses
            fprintf('     data state %1.0f(%1.0f): %1.0fx%1.0f, ',i,j,datsiz(1),datsiz(2))
        end
        switch lam
            case -3
                fprintf('healthy/immobile/dead data.')
            case -2
                fprintf('survival data (dose-response).')
            case -1
                fprintf('survival data.')
            case 0
                if prod(datsiz) == 0
                    fprintf('no data.')
                else
                    fprintf('continuous data, log transf.')
                end
            case 1
                fprintf('continuous data, no transf.')
            otherwise
                fprintf('continuous data, power %1.2f transf.',lam)
        end
        if isfield(glo,'wts') % weights have been applied
            if n_D > 1
                fprintf(' Weight: %1.4g.',glo.wts(j,i))
            else
                fprintf(' Weight: %1.4g.',glo.wts(i))
            end
        end
        if isfield(glo,'var') % variances have been supplied
            if n_D > 1
                fprintf(' Var: %1.4g.',glo.var(j,i))
            else
                fprintf(' Var: %1.4g.',glo.var(i))
            end
        end
        fprintf('\n')
    end
end

switch opttype
    case 1
        fprintf('   Search method: Nelder-Mead simplex direct search, %g rounds. \n',simno)
    case 2
        fprintf('   Search method: Simulated annealing plus simplex search, stop temp. %g. \n',anntemp)
    case 3
        fprintf('   Search method: Swarm optimisation plus simplex search, %g particles, %g iter. \n',swno,swit)
    case 4
        fprintf('   Search method: Parameter-space explorer (see CIs above). \n')
    case 5
        fprintf('   Search method: alternative Nelder-Mead simplex direct search (nelmin), %g rounds. \n',simno)
end

if EXITFLAG == 1
    disp('     The optimisation routine has converged to a solution')
else
    disp('     No convergence, optimisation stopped on maximum number of iterations')
end
if ~isinf(out_iter)
    fprintf('     Total %1.0f simplex iterations used to optimise. \n',out_iter)
end
fprintf('     Minus log-likelihood has reached the value %#1.2f (AIC=%#1.2f). \n',FVAL,2*length(pfit)+2*FVAL)

disp('=================================================================================')
nfields = length(names);
for i = 1:nfields % display all parameters on screen
    if pmat(i,2) == 0 % parameter is not fitted
        fprintf('%-6s %10.4g (fit: %1.0f, initial: %2.4g) ',names{i},p(i),pmat(i,2),p0(i))
    else % parameter is fitted: also print extra zeros on estimate
        fprintf('%-6s %#10.4g (fit: %1.0f, initial: %2.4g) ',names{i},p(i),pmat(i,2),p0(i))
    end
    if ~isempty(glo2.pri)
        pri = glo2.pri;
        if isfield(pri,names{i}) % is there a prior for this one?
            a = eval(['pri.',names{i}]); % extract prior info
            switch a(1)
                case 1
                    fprintf('\t prior: beta (%g,%g)',a(2),a(3))
                case 2
                    fprintf('\t prior: tri  (%g-%g-%g)',a(2),a(4),a(3))
                case 3
                    fprintf('\t prior: norm (mn %g,sd %g)',a(2),a(3))
                case 4
                    fprintf('\t prior: logn (lnmn %g,lnsd %g)',a(2),a(3))
                otherwise
                    fprintf('\t prior: unknown distribution ...')
            end
        else
            fprintf('\t prior: unif (%g-%g)',pmat(i,3),pmat(i,4))
        end
    end
    fprintf('\n')
end

ind_log_fit = find(pmat(:,5)==0 & pmat(:,2)==1); % find fitted parameters that are on log scale
if ~isempty(ind_log_fit)
    disp('=================================================================================')
    
    if length(ind_log_fit)==1
        fprintf('Parameter ')
    else
        fprintf('Parameters ')
    end
    for i = 1:length(ind_log_fit)
        fprintf('%s',names{ind_log_fit(i)})
        if i ~= length(ind_log_fit)
            if i == length(ind_log_fit)-1
                fprintf(' and ')
            else
                fprintf(', ')
            end
        end
    end
    if length(ind_log_fit)==1
        fprintf(' is fitted on log-scale.\n')
    else
        fprintf(' are fitted on log-scale.\n')
    end
end

disp('=================================================================================')
if ~isempty(zvd)
    disp('Final values for zero-variate data points included:')
    for i = 1:length(namesz) % go through contents of the global zvd
        a = eval(['zvd.',namesz{i}]); % extract zero-variate info
        fprintf('%-6s %#10.4g (data point: %1.4g, sd: %1.4g) \n',namesz{i},a(3),a(1),a(2))
    end
    disp('=================================================================================')
end
% for analysis of toxicity data, we can print more info for the fit
if isfield(glo,'moa')
    disp(['Mode of action used: ',num2str(glo.moa)])
    if isfield(glo,'feedb')
        disp(['Feedbacks used     : ',num2str(glo.feedb)])
    end
    disp('=================================================================================')
end

disp(['Time required: ' secs2hms(toc)])
diary off
