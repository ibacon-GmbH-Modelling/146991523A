function [data_out,S0,w_out] = recalc_data(data_in,w_in,lam)

% Usage: [data_out,S0,w_out] = recalc_data(data_in,w_in,lam)
%
% Function to recalculate the data for a SINGLE treatment for plotting. If
% entering replicates, it will combine replicates into a single value with
% CI. The function also deals with NaNs, weights/outliers (sub-lethal
% endpoints) and missing/removed animals (survival).
% 
% For survival, the CI is the Wilson score interval (also used in
% openGUTS). For sub-lethal data, weights are used to calculate a weighted
% mean. The CI on the weighted mean is approximated by two times the SE of
% the weighted mean, as calculated according to Madansky and Alexander
% (following WinCross/Quantum). See:
% http://www.analyticalgroup.com/download/WEIGHTED_MEAN.pdf. For this to
% work properly, weights need to represent the number of observations on
% which a data point is based (i.e., the number of animals in the replicate
% on which a mean is based).
%
% This function is called in plot_tktd and in calc_and_plot. Large parts of
% the code is based on openGUTS prepare_data.
%
% Author     : Tjalling Jager
% Date       : March 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM.

if lam == -1 || lam == -3 % than it is a survival data set
    
    indnan = isnan(data_in); % remember where the NaNs in the data are
    
    % below, the vector/matrix with missing animals is used (cumulative)
    % for this treatment only. Note that first row in w_ic must be zeros
    S0 = ones(size(data_in,1),1) * data_in(1,:) - cumsum(cat(1,w_in(1,:),w_in(1:end-1,:)),1); % initial number of animals over time, corrected for missing animals
    S0(indnan) = 0; % where the data is NaN, make the number of animals zero
    S0 = sum(S0,2); % total number of animals over all replicates, corrected for NaNs and missing animals
    
    s_ic = sum(data_in,2,'omitnan'); % for survival data, we'll simply sum the survivors, ignoring the NaNs!
    % Without missing animals, this can be used to calculate survival
    % probabilities directly, but even with missing animals this is needed
    % to calculate deaths.
    
    if sum(w_in(:)) == 0 % if there are NO missing animals, things are simple!
        w_in  = sum(w_in,2); % rather superfluous as they are zero anyway
        ps_ic = s_ic ./ S0 ; % vector with survival probabilities
        % Note: where all replicates are NaN, S0 will be zero, so the
        % probability to survive will become NaN again. The Wilson score
        % will also become NaN. When only some replicates are NaN, this
        % works fine: the probability is calculated based on the replicates
        % that ARE there.
    else
        % For missing/removed animals, the calculation is a bit more
        % cumbersome. I initially thought that it was okay to use S0 to
        % calculate the probabilities: if you remove animals, that simply
        % means you started out with less. However, it does not work that
        % way. Perhaps easiest to illustrate with an example. Suppose we
        % have survival data [100;50;10] with 40 animals missing after the
        % second observation. Initial nr of animals is thus for the two
        % intervals 100 and 60. The survival probability would then, with
        % S0, be [1;0.5;0.17], which is incorrect as there is no mortality
        % at all in the second interval. The correct way was already used
        % in DEBtoxM and in BYOM, so I just needed to rediscover it. You
        % need to calculate *conditional* survival probabilities for each
        % interval and multiply them. The implementation is different from
        % previous versions to make it more tractable (and as we keep the
        % NaNs in the data set).
        for i = 1:size(data_in,1) % run through time points
            if any(indnan(i,:)) % if there is at least one NaN on this time point
                data_in(i,indnan(i,:)) = data_in(i-1,indnan(i,:))-w_in(i-1,indnan(i,:));
            end
            % replace NaN in data set by previous observation, minus
            % removed animals there; this assumes that no animals die in
            % the NaN intervals, which is okay as we remove these intervals
            % later anyway. This also works when there are more NaNs in the
            % same column, at consecutive time points (that is why I think
            % this needs a for loop). This is a bit awkward, but is
            % difficult to avoid as there may be a NaN in one or more (or
            % all) replicates at some time point.
        end
        
        Ninit = data_in - w_in; % initial number of animals at start of each interval
        Nend  = [data_in(2:end,:);zeros(1,size(data_in,2))]; % living animals at the end of each interval
        % Sum the replicates! I keep Ninit and Nend as well, as they may be
        % useful in the future to work with NaNs in some way.
        Ninit_tot = sum(Ninit,2);
        Nend_tot  = sum(Nend,2);
        
        Pcond = [1 ; Nend_tot(1:end-1,:) ./ Ninit_tot(1:end-1,:)]; % conditional probability to survive each interval
        Pcond(Ninit_tot==0) = 0; % where initial number is zero, probability is zero (avoid NaNs)
        ps_ic = cumprod(Pcond); % unconditional survival probability as product across time
        ps_ic(sum(indnan,2)~=0) = NaN; % place the NaNs back where they belong
        % This combines the replicates, but sets a NaN also where ONE of
        % the replicates gives a NaN ... this might be best as there seems
        % to be no good way to calculate a survival probability when one of
        % the replicates is missing.
        
    end
    
    % Wilson score interval on data probabilities.
    % https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
    z = 1.96; % critical value of normal distribution for 95% CI
    n = S0; % number of animals to start with
    a = (ps_ic + z^2./(2*n))./(1+z^2./n); % middle of confidence interval (element-wise calculation on the vector <ps_ic>)
    b = z./(1+z^2./n) .* sqrt(ps_ic.*(1-ps_ic)./n + z^2./(4*n.^2)); % range of interval (element-wise calculation on the vector <ps_ic>)
    b(1) = 0; % at t=0, there is no uncertainty, so remove the range
    a(1) = 1; % and make middle of range 1
    % Note: the interval itself (the min and max) are given by a-b and
    % a+b, which is done in the next block.
    
    % I now leave the NaNs in there!
    indnan_ps = isnan(ps_ic);
    
    % Can we calculate deaths when combining replicates? Only
    % if there are no NaN, or both replicates have the NaNs at
    % the same time point ... I check this on the means of the
    % logicals returned by isnan:
    chck_nans = mean(indnan,2); % this is 0 when there are no NaNs, and 1 when there are NaNs in all replicates
    if all(chck_nans == 1 | chck_nans == 0) % if they are all 0 or 1, we can calculate deaths as usual
        % Problem is that someone may want to remove animals at a time when
        % there is no observation on survival numbers ... This should
        % already generate an error in prelim_checks.
        if sum(w_in(indnan_ps)) ~= 0
            error('There are missing/removed animals at a time point where survival is NaN ... this needs more work!')
        end
        d_ic = nan(size(s_ic));
        d_ic(~indnan_ps) = -diff([s_ic(~indnan_ps);0]) - w_in(~indnan_ps);
        % this places the deaths at the same time vector, with NaN for deaths
        % at a time point where there is a NaN for survival
    else % then there is at least one replicate with a NaN(s) and non-NaN(s) at same time point
        d_ic = nan(size(s_ic)); % need to figure out how to calculate deaths ... 
        % This may not be possible as a NaN in one replicate will make it
        % impossible to specify the interval for the deaths, but should
        % perhaps be possible to do it crudely.
        %
        % Note: a warning will be produced on screen by plot_tktd and
        % calc_and_plot when attempting to plot means with this data set.
    end
    
    % Combine all vectors into the overall object, as a nice matrix:
    % surv. prob., CI min, CI max, survivors, nr. deaths
    data_out = [ps_ic max(0,a-b) min(1,a+b) s_ic d_ic]; % creates a matrix
    data_out(indnan_ps,:) = NaN; % make sure all entries are NaN where observation is NaN
    w_out = w_in;
    
else % it is a continuous data set
    
    % Note: below, the mean is calculated in a somewhat odd manner (without
    % using Matlab's mean function). This is done to accommodate the fact
    % that weights may differ between replicates; we have to do a weighted
    % mean.
    %
    % Note: we want to keep outliers (don't make their mean a NaN), and we
    % need to account for fact that weights and NaNs might differ between
    % replicates. If all replicates are outliers, the mean is calculated
    % using all points (but the mean will be recognised as outlier by the
    % plotting routines). If some replicates are outliers, but not all of
    % them, the mean is calculated from the values that are NOT outliers,
    % and the point will NOT be marked as outlier in the plotting routines.
    
    w_in(isnan(data_in)) = NaN; % make sure that weights are NaN where data is NaN (this is needed to make the calculations work)
    ind_zero    = find(sum(w_in,2,'omitnan') == 0); % indices to outliers (across all replicates), that we want to treat separately
    ind_nonzero = find(sum(w_in,2,'omitnan') ~= 0); % indices to rows across which to calculate a mean
    
    if lam == 0 % than do log-transformation before calculating mean and SE
        data_in(~isnan(data_in)) = max(data_in(~isnan(data_in)),1e-10); % make sure that there are no zeros in there
        data_in = log(data_in); % take the LN of the data
    else % than do a power transformation (none if lam=1) before taking residuals
        data_in = data_in.^lam; % transform data
    end
    
    data_mn = nan(size(data_in,1),1); % initiate a vector to catch the means
    data_mn(ind_zero)    = mean(data_in(ind_zero,:),2,'omitnan'); % for the outliers, calculate a mean in the normal way (there are no weights anyway)
    data_mn(ind_nonzero) = sum(data_in(ind_nonzero,:) .* w_in(ind_nonzero,:),2,'omitnan') ./ sum(w_in(ind_nonzero,:),2,'omitnan');
    % for the regular points, calculate a weighted mean
    
    data_in(w_in == 0)   = NaN; % make outliers NaN in the data; this way, outliers don't attribute to CI!
    data_in(isnan(w_in)) = NaN; % make sure that data are NaN where weights is NaN (NaN weights also don't contribute)
    % this also means that outliers with ind_zero do not get a CI ...
    
%     % this is the regular, unweighted, standard error
%     data_se = std(data_in,0,2,'omitnan')./sqrt(sum(~isnan(data_in),2)); % calculate standard error and ignore nans
    
    % This is the SE of the weighted mean, calculated according to Madansky
    % and Alexander, following the WinCross/Quantum approach.
    % See: http://www.analyticalgroup.com/download/WEIGHTED_MEAN.pdf. This
    % also works when there are no weights, but then all w_in should be 1!
    data_sd  = std(data_in,0,2,'omitnan'); % standard deviation from the data
    data_var = data_sd .^2; % variance of the data
    b        = ((sum(w_in,2,'omitnan')).^2) ./ sum(w_in.^2,2,'omitnan'); % "effective base"
    var_mean = data_var ./ b; % variance of the weighted mean (is SE^2)
    data_se  = sqrt(var_mean); % standard deviation of the weighted mean (is SE)
    data_ci  = [data_mn-2*data_se data_mn+2*data_se]; % CI on transformed scale (2 times SE)
        
    % back-transform mean and CI when needed
    if lam == 0 % then we did log-transformation before calculating mean and SE
        data_mn = exp(data_mn);
        data_ci = exp(data_ci);
    else % then we did power transformation before calculating mean and SE
        data_mn = data_mn.^(1/lam);
        data_ci = data_ci.^(1/lam);
    end
    
    data_out = [data_mn data_ci]; % collect output in matrix
    w_out    = sum(w_in,2,'omitnan'); % sum the weights, take NaNs as zero weight
    S0       = [];
    
end
