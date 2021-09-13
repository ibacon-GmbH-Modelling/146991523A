function [minloglik,GoodFit,zvd] = transfer(pfit,pmat,WRAP)

% Usage: [minloglik,GoodFit,zvd] = transfer(pfit,pmat,WRAP)
%
% This is a transfer function which calculates the minus log-likelihood
% function for a certain set of parameters.
%
% This transfer function allows certain parameters to be fixed, while
% others are fitted. The free parameters are in <pfit>, all parameters in
% <p_mat>, which contains in the second column a 1 where the parameter is
% to be fitted. <X0mat> is the matrix with initial states.
%
% This function is used with <fminsearch> (and the other optimisation and
% sampling routines), and when calculating CRs and CIs (such as in
% <calc_proflik.m>).
%
% Author     : Tjalling Jager 
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

% Unpack the structure WRAP into its components (previously global)
DATA  = WRAP.DATA;
W     = WRAP.W;
DATAx = WRAP.DATAx;
Wx    = WRAP.Wx;
X0mat = WRAP.X0mat;
glo   = WRAP.glo;
glo2  = WRAP.glo2;

sameres = 0; % use one overall nuisance residual sd for all datasets for one state (this needs further testing!)
Fsdmin  = 0; % option to force a minumum value of the SSQ to avoid focus on perfect fit (try 0.001)

% extract parameters from the general globals, glo and glo2
n_D  = glo2.n_D;
n_X  = glo2.n_X;
n_X2 = glo2.n_X2;
ttot = glo2.ttot;
ctot = glo2.ctot;

p = pmat(:,1);           % copy the whole parameter vector to p
p(pmat(:,2)==1) = pfit ; % replace values in p with the values to be fit (pfit)

% NOTE: code below is modified to check whether the parameter is within
% bounds on the log-transformed value. Reason is that 10^(log10(x))-x is
% not equal to x, but 3.5527e-15. This means that when an element of p_tmp
% would be exactly on the bound, it will still be *seen* as outside.

% put parameter bounds for parameters fitted on log scale also on log scale
pmat(pmat(:,5)==0,[3 4]) = log10(pmat(pmat(:,5)==0,[3 4]));
% if any parameter is outside the bounds, return +inf as criterion
if any(p<pmat(:,3) | p>pmat(:,4))
    minloglik = +inf; % than make the minloglik infinite (very bad fit!)
    return; % and return to the function calling this one
end
% put parameters that need to be fitted on log scale back on normal scale
p(pmat(:,5)==0) = 10.^(p(pmat(:,5)==0));
% pmat needs to be on normal scale again for call_deri

ndata   = n_X * n_D; % total number of data sets
GoodFit = nan(ndata,1); % initialise Goodfit with NaNs

% initialise the matrices that will collect the model output
Xcoll = cell(n_X,1); % create empty cell array
for i = 1:n_X     % loop over the data sets
    Xcoll{i} = zeros(length(ttot),length(ctot)); % and fill an Xcoll with zeros
end
if n_X2 > 0 % if there are additional univariate data ...
    Xcoll2x = cell(n_X2,1); % create empty cell array
    Xcoll2y = cell(n_X2,1); % create empty cell array
%     for j = 1:n_X2 % now loop over the additional data sets
%         Xcoll2x{j} = zeros(length(ttot),length(ctot)); % and fill an Xcoll2x with zeros
%         Xcoll2y{j} = zeros(length(ttot),length(ctot)); % and fill an Xcoll2y with zeros
%     end
% % commented out as extra data may not have time on x-axis
end

par = packunpack(2,0,p,WRAP); % transform parameter matrix into a structure

%% Calculate model results
% =========================================================================

for i = 1:length(ctot) % run through all of our concentrations
    % initial states to start the solver with, from X0mat
    X0 = X0mat(2:end,X0mat(1,:) == ctot(i)); % read initial states from X0mat using "concentration" vector
    
    % Note: glo is here handed over with the call to call_deri. That means
    % that it does not have to be global in call_deri. 
    
    [Xout,~,Xout2,zvd] = call_deri(ttot,par,[ctot(i);X0],glo); % use call_deri.m to provide the output per concentration
    % Note: zvd is glo.zvd with an added (third) element: the estimated
    % value based on the parameters in par (and possibly other things).
   
    if size(Xout,1) == length(ttot)
        for j = 1:n_X % now loop over the state variables
            Xcoll{j}(:,i) = Xout(:,j); % and collect the results in the right part of Xcoll
        end
    else % call_deri did not return all required time points, so something went wrong in calculating the derivatives
        minloglik = +inf; % than make the minloglik infinite (very bad fit!)
        % error('This should not happen anymore ...')
        return; % and return to the function calling this one
    end
    for j = 1:n_X2 % now loop over the additional data sets
        Xcoll2x{j}(:,i) = Xout2{j}(:,1); % collect the new x-values
        Xcoll2y{j}(:,i) = Xout2{j}(:,2); % collect the new y-values
    end
end

%% Calculate likelihood
% =========================================================================

loglik = zeros(ndata,1); % initialise a vector for the log-likelihoods with zeros
find_data = reshape(1:ndata,n_D,n_X); % helper to find which state a data set is
if n_D > 1 && sameres == 1
    i_rem = 1; % perhaps also initialise rem_ssq?
end

for i = 1:ndata % loop over the data sets
    
    if size(DATA{i},1) == 1 % if the DATA is a row vector (only scenarios, or just a zero), ignore it
        loglik(i) = 0;
    else
        
        D   = DATA{i}(2:end,2:end); % extract the correct data matrix
        tD  = DATA{i}(2:end,1);     % time vector in this data set
        cD  = DATA{i}(1,2:end);     % concentration vector in this data set
        lam = DATA{i}(1,1);         % -1 for survival data, 0-1 for continuous data (data transformation when lam<1)
        w   = W{i};                 % weight factors are by default set to 1 in start_calc
        
        [~,locC]=ismember(cD,ctot); % locC is location in ctot for common values with cD
        if any(locC==0) % it is possible that the model scenarios do not cover all scenarios in the data
            D(:,locC==0)  = []; % remove columns of data that we do not want to fit (they are not in ctot)
            w(:,locC==0)  = []; % remove columns of weights that we do not want to fit (they are not in ctot)
            locC(locC==0) = []; % also remove the zeros from the location vector
        end
        [~,locT]=ismember(tD,ttot); % locT is location in ttot for common values with tD
        if any(locT==0) % if we set glo.Tinit, it is possible that the model times do not cover all scenarios in the data
            D(locT==0,:)  = []; % remove columns of data that we do not want to fit (they are not in ctot)
            w(locT==0,:)  = []; % remove columns of weights that we do not want to fit (they are not in ctot)
            locT(locT==0) = []; % also remove the zeros from the location vector
        end
        
        if ~isempty(D) % D will be made empty when the scenarios do not match the data at all
            
            if n_D == 1
                M = Xcoll{i}(locT,locC); % model values at the concentration and timepoint of the data
                % this also works when there are replicated data; the model values
                % are simply copied to the correct columns!!
                i_st = i; % then i is also the number of the state we currently process
            else
                [~,i_st] = find(find_data==i); % find the state that data set i belongs to
                M = Xcoll{i_st}(locT,locC); % model values at the concentration and timepoint of the data
            end
            
            % Binomal distribution for dose-response curve fit (package doseresp)
            if lam == -2 % we have survival data in dose-response context
                % w now carries the starting animals
                loglik(i) = 0;
                for i_d = 1:numel(D) % run through elements of D, and calculate binomial probability
                    if ~isnan(D(i_d)) % simply ignore NaNs
                        % make sure survival probability is not zero or one (this gives trouble taking log)
                        p = min(1-1e-10,M(i_d));
                        p = max(1e-10,p);
                        logprob = D(i_d)*log(p) + (w(i_d)-D(i_d))*log(1-p); % log-lik up to a constant
                        loglik(i) = loglik(i) + logprob; % add this treatment to the previous ones
                    end
                end
                
                % Do silly r-square for survival.
                S = D./w; % survival probability from the data set
                ind_ok = ~isnan(S); % are there any NaNs?
                S = S(ind_ok); % remove them from the data probabilities
                M = M(ind_ok); % and also from the model probabilities
                res = (S - M) ; % residuals for survival probability
                res_tot = (S) - mean(S) ; % residuals from mean of data
                GoodFit(i) = 1 - (res' * res)/(res_tot' * res_tot); % R-square, NOT adjusted for no. of pars.
                
            end
            
            if lam == -3 % we have multi-state quantal data in unconditional MN setting
                
                % make sure probability is not zero (this gives trouble taking log)
                p = max(1e-10,M);
                
                if sum(w(:)) ~= 0
                    error('Cannot deal yet with missing/removed animals for immobility!')
                end
                S0 = ones(size(D,1),1) * D(1,:); % initial number of animals over time
                
                % % Attempt to work with missing animals; this needs more thought
                % S0 = ones(size(D,1),1) * D(1,:) - cumsum(cat(1,w(1,:),w(1:end-1,:)),1); 
                % % initial number of animals over time, corrected for
                % % missing animals (do this before removing row 1, as row 1
                % % of D contains the initial number of individuals!)
                
                % we need to correct the data for dead, immobile and
                % immobile+dead; the first row is the total number of
                % individuals to start with (to allow calculating a
                % fraction for plotting), but for the likelihood we need to
                % make them zeros
                if i_st ~= glo.loc_h
                    D(1,:) = 0;
                end
                
%                 % Likelihood contributions, but ignore response at t=0?
%                 lik_contr = D .* log(p); % likelihood contribution for each observation
%                 lik_contr = lik_contr(:); % make it into a column vector
%                 loglik(i) = sum(lik_contr,1,'omitnan'); % sum all contributions, ignore NaNs
                
                % some more work is needed to calculate r-square ...                
                S = D ./ S0; % probabilities from the data, corrected for missing animals
                S = S(:);
                M = M(:);
                M(isnan(S)) = []; % remove places where NaNs in data set
                S(isnan(S)) = []; % remove places where NaNs in data set
                
                res        = (S - M) ; % residuals for survival probability
                res_tot    = (S) - mean(S) ; % residuals from mean of data
                GoodFit(i) = 1 - (res' * res)/(res_tot' * res_tot); % R-square, NOT adjusted for no. of pars.
                
                % For now, use SSQ instead of multinomial likelihood! The
                % reason is that some data sets show just a few immobiised
                % individuals after substantial recovery time, which
                % heavily affects the likelihood function. This needs more
                % thought, but for now the SSQ will do.
                loglik(i) = -(length(S)/2) * log(res' * res); % loglik from sum of squares
                
            end
            
            if lam == -1 % than we have survival data, in multinomial context
                
                if sum(w(:)) == 0
                    % To calculate a r-square (which is a bit silly for
                    % survival data, though it can be called a 'model
                    % efficiency'), we need the survival fraction from the
                    % data. This works for NaNs but not for missing data
                    % (yet; this could be done as in recalc_data but that
                    % would be calculation intensive here as transfer is
                    % called many times). However, plot_tktd will provide
                    % an r-square when possible.
                    % 
                    % We can also remove the r-square calculation here to
                    % win some time ...
                    S0 = ones(size(D,1),1)*D(1,:); % initial number of animals over time, corrected for missing animals
                    S  = D ./ S0; % survival probability from data
                    S  = S(:); % turn it into a row vector
                    indnan = isnan(S); % location of NaNs
                    S(indnan) = []; % remove the NaNs
                    M_tmp = M(:); % turn it into a column vector
                    M_tmp(indnan) = []; % remove the NaNs where S is NaN
                    % Note that the r-square calculation includes the
                    % response at t=0, which is fitted perfectly anyway. I
                    % think this is better than excluding it.
                    res        = (S - M_tmp) ; % residuals for survival probability
                    res_tot    = (S) - mean(S) ; % residuals from mean of data
                    GoodFit(i) = 1 - (res' * res)/(res_tot' * res_tot); % R-square, NOT adjusted for no. of pars.
                end
                
                % The following part works also with NaNs and missing data.
                % The trick is to look at one treatment at a time, and
                % remove the positions with NaNs. This does not work when
                % there are missing data entered at a time point where
                % observed survivors is NaN, but this yields an error in
                % prelim_checks already. 
                % 
                % Note: the implementation has changed now in v. 5.0. The
                % new version is more readable, and equally fast (perhaps
                % slightly slower for GUTS analyses without missing animals
                % or NaNs).
                    
                logliktmp = 0;
                indnan    = isnan(D); % indices to entries with NaNs
                for i_t = 1:size(D,2) % run through treatments in data set
                    indnan_i      = indnan(:,i_t); % NaN indices for this treatment
                    D_i           = D(:,i_t);   % data for this treatment
                    D_i(indnan_i) = [];         % remove the entries with NaNs
                    M_i           = M(:,i_t);   % model points for this treatment
                    M_i(indnan_i) = [];         % remove the entries with NaNs
                    w_i           = w(:,i_t);   % missing animals for this treatment
                    w_i(indnan_i) = [];         % remove the entries with NaNs
                    
                    Ndeaths = -diff([D_i;0])-w_i; % numbers of deaths, corrected for missing animals
                    Mdeaths = -diff([M_i;0]);     % conditional probabilities of deaths
                    Mdeaths = max(Mdeaths,1e-50); % otherwise we get problems taking the logarithm
                    M_i     = max(M_i,1e-50);     % otherwise we get problems taking the logarithm
                    if any(Ndeaths<0)
                        error('Negative deaths discovered. Please check the weights matrix in the data set.')                        
                    end
                    
                    logliktmp = logliktmp + Ndeaths'*log(Mdeaths) + w_i'*log(M_i);  % and calculate the log likelihood
                    % note that the removed or missing animals contribute too by
                    % the information that they WERE alive up till some time point
                end
                loglik(i) = logliktmp;
                 
            end
                        
            if lam >= 0 % than we have continuous response data
                
                D = D(:); % put data in a single vector
                M = max(0,M(:)); % and the same for the model results
                w = w(:); % and the weight factors
                ind_fin = (isfinite(D)==1 & w~=0); % find where data is a nice finite number (so not NaN or zero)
                % note: don't use points with zero weight in the index
                
                % if sum(ind_fin) == 0 % there is no useful data at all in this data set
                %    loglik(i) = 0;
                %    GoodFit(i) = NaN;
                %    break % and move out of the for loop
                % end                    
                
                % remove the ones that are NaN or weight zero (30/1/2020
                % cleaned the code a bit by removing elements from D, M and
                % w, rather than using the ind_fin in all calculations).
                D = D(ind_fin);
                M = M(ind_fin);
                w = w(ind_fin);
                
                if lam == 0 % than do log-transformation before taking residuals
                    D   = max(D,1e-10);     % zeros lead to problems ... this is nasty when the data elements are tiny ...
                    M   = max(M,1e-10);
                    res = log(D) - log(M) ; % residuals from transformed data
                    mn  = mean(log(D));     % mean value of data after transformation
                    res_tot = log(D) - mn ; % residuals from mean of data
                else % than do a power transformation (none if lam=1) before taking residuals
                    res = D.^lam - M.^lam ; % residuals from transformed data
                    mn  = mean(D.^lam);     % mean value of data after transformation
                    res_tot = D.^lam - mn ; % residuals from mean of data
                end
                wssq  = res' * diag(w) * res ; % weigthed sum of squares
                wssq2 = res' * (diag(w)).^2 * res ; % weigthed sum of squares (weighted with n squared)
                n     = sum(ind_fin); % number of data points (already excluding points with zero weight)
                N     = sum(w);       % total nr of individuals used
                wssq_tot   = res_tot' * diag(w) * res_tot ; % weigthed sum of squares
                GoodFit(i) = 1 - wssq/wssq_tot; % R-square, NOT adjusted for no. of pars.
                % GoodFit(i) = 1 - (wssq/(n-length(pfit)-1))/(wssq_tot/(n-1); % R-square, adjusted for no. of pars. 

                % From Wikipedia: Adjusted R2 does not have the same
                % interpretation as R2�while R2 is a measure of fit,
                % adjusted R2 is instead a comparative measure of
                % suitability of alternative nested sets of explanators.
                
                if wssq_tot == 0 % there is no variation in the data ...
                    GoodFit(i) = NaN; % make it a NaN (otherwise it is -INF)
                end
                
                if ~isfield(glo,'var')
                    % set a minimum for the wsssq, so that sd of residuals is >
                    % a factor Fsdmin (top of function) of the average of the data
                    % this should avoid focus on perfect fits
                    wssq2 = max(wssq2,N*(Fsdmin * mn)^2);
                    wssq  = max(wssq,n*(Fsdmin * mn)^2);
                    % Note: watch out when fitting means with variable number
                    % of individuals per mean (not exactly correct, so
                    % discontinuity may occur using max).
                    
                    if n_D == 1 || sameres == 0 % for one data set per state, simply calculate log likelihod
                        % calculate log-likelihood for replicated data with sd as nuisance parameter
                        loglik(i) = -(n/2) * log(wssq2) - N*wssq/(2*wssq2);
                    else % if there are more data sets that we need to combine ...
                        rem_ssq(i_rem,:) = [i_st i wssq wssq2 n N];
                        i_rem = i_rem + 1;
                    end
                else
                    % calculate log-likelihood for replicated data with sd provided in global
                    % the variance must be for individuals, and after transformation!
                    % If variance of mean is known, enter variance times number
                    % of replicates. This makes sure that differences between
                    % number of replicates can be dealt with properly.
                    loglik(i) = -1/(2*glo.var(i)) * wssq;
                end
            end
        end
    end
end

if n_D > 1 && sameres == 1 % in case we want to have a single residual sd per state ...
    rem_un = unique(rem_ssq(:,1));
    for i = 1:length(rem_un) % go through unique states in the remembered matrix with ssq's and n's
        ind_i = find(rem_ssq(:,1)==rem_un(i));
        rem_i = sum(rem_ssq(ind_i,3:end),1);
        % calculate a log-likelihood from the summed results for this state
        % and place it in loglik at first position for this state
        loglik(rem_ssq(ind_i(1),2)) = -(rem_i(3)/2) * log(rem_i(2)) - rem_i(4)*rem_i(1)/(2*rem_i(2));
    end
end

if isfield(glo,'wts') % if weights for the data sets are specified ...
    loglik = loglik .* glo.wts(:); % use weights for each data type
    % Note: when using a common residual error for multple data sets within
    % a state, only the first weight for that state will be used!
    loglik(isnan(loglik) & glo.wts(:)==0) = 0; 
    % This is to make sure that a state can be NaN as long as it has weight zero!
end

logpriprob = 0; % initiate prior probability
if ~isempty(glo2.pri) % note that par is on normal scale, so priors are defined at normal scale
    logpriprob = calc_prior(par,glo2); % call function to calculate prior prob.
end

logzvdprob = 0; % initiate zero-variate data probability
if ~isempty(zvd)
    logzvdprob = calc_zerovar(zvd,glo2); % call function to calculate zero-variate data prob.
end

% Next, calculate the log-likelihood for any additional data sets
loglik2 = zeros(n_X2,1); % initialise a vector for the log-likelihoods with zeros
for i = 1:n_X2 % loop over the additional data sets (not done when n_X2 = 0)
    
    if size(DATAx{i},1) == 1 % if the DATA is a row vector (only scenarios, or just a zero), ignore it
        loglik2(i) = 0;
    else
        
        D   = DATAx{i}(2:end,2:end); % extract the correct data matrix
        xD  = DATAx{i}(2:end,1);     % x vector in this data set
        cD  = DATAx{i}(1,2:end);     % concentration vector in this data set
        lam = DATAx{i}(1,1);         % -1 for survival data, 0-1 for continuous data (data transformation when lam<1)
        w   = Wx{i};                 % weight factors are by default set to 1 in start_calc

        [~,locC]=ismember(cD,ctot); % locC is location in ctot for common values with cD
        if any(locC==0) % it is possible that the model scenarios do not cover all scenarios in the data
            D(:,locC==0)  = []; % remove columns of data that we do not want to fit (they are not in ctot)
            w(:,locC==0)  = []; % remove columns of weights that we do not want to fit (they are not in ctot)
            locC(locC==0) = []; % also remove the zeros from the location vector
        end
        
        if ~isempty(D) % D will be made empty when the scenarios do not match the data at all
            
            Mx = Xcoll2x{i}(:,locC);    % model x-values at the concentration of the data
            My = Xcoll2y{i}(:,locC);    % model y-values at the concentration of the data
            % this also works when there are replicated data; the model values
            % are simply copied to the correct columns!!
            M = nan(length(xD),size(Mx,2)); % initialise with NaNs
            for j = 1:size(Mx,2) % run through scenarios
                M(:,j) = interp1(Mx(:,j),My(:,j),xD,'linear','extrap'); % linear inter- and extrapolation
                % to get model values at the exact x-values of the data set
            end
            
            D = D(:); % put data in a single vector
            M = max(0,M(:)); % and the same for the model results
            w = w(:); % and the weight factors
            ind_fin = (isfinite(D)==1 & w~=0); % find where data is a nice finite number (so not NaN or zero)
            % note: don't use points with zero weight in the index
            
            % remove the ones that are NaN or weight zero (30/1/2020
            % cleaned the code a bit by removing elements from D, M and w,
            % rather than using the ind_fin in all calculations).
            D = D(ind_fin);
            M = M(ind_fin);
            w = w(ind_fin);
                
            if lam == 0 % than do log-transformation before taking residuals
                res = log(D+1e-10) - log(M+1e-10) ; % residuals from transformed data
                mn  = mean(log(D+1e-10)); % mean value of data after transformation
            else % than do a power transformation (none if lam=1) before taking residuals
                res = D.^lam - M.^lam ; % residuals from transformed data
                mn  = mean(D.^lam); % mean value of data after transformation
            end
            wssq  = res' * diag(w) * res ; % weigthed sum of squares
            wssq2 = res' * (diag(w)).^2 * res ; % weigthed sum of squares (weighted with n squared)
            n = sum(ind_fin); % number of data points
            N = sum(w);       % total nr of individuals used
        
            % set a minimum for the wsssq, so that sd of residuals is >
            % a factor Fsdmin (top of function) of the average of the data
            % this should avoid focus on perfect fits
            wssq2 = max(wssq2,N*(Fsdmin * mn)^2);
            wssq  = max(wssq,n*(Fsdmin * mn)^2);
            % Note: watch out when fitting means with variable number
            % of individuals per mean (not exactly correct, so
            % discontinuity may occur using max).
            
            % calculate log-likelihood for replicated data with sd as nuisance parameter
            loglik2(i) = -(n/2) * log(wssq2) - N*wssq/(2*wssq2);
           
        end
    end
end

if isfield(glo,'wts2') % if weights for the extra data sets are specified ...
    loglik2 = loglik2 .* glo.wts2(:); % use weights for each data type
end

% Calculate the overall log likelihood
minloglik = -1 * (sum(loglik)+logpriprob+logzvdprob+sum(loglik2)) ; % makes it negative. This is needed as fminsearch is looking for a MINIMUM,
% and we need to find the MAXIMUM of the likelihood function.

if isinf(minloglik) % the minloglik should not be calculated as infinite ...
    error('The min-log-likelihood is calculated to be infinite ... check data set and derivatives or debug from transfer.m')
end

if isnan(minloglik) || ~isreal(minloglik) % when the loglik calculation returns NaN or complex numbers ...
    minloglik = +inf; % give it a really bad likelihood value so the optimisation ignores it
end

