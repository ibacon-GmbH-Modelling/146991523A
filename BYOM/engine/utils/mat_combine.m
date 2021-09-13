function Ctot = mat_combine(estim,varargin)

% Usage: Ctot = mat_combine(estim,varargin)
%
% A smart little tool to combine data sets with different time vectors into
% one data matrix. The option <estim> can be used to fill missing spaces
% with NaN (estim=0) or using splining to interpolate them (estim=1).
% As input, the function needs matrices with concentrations in the first
% row and time in the first column. The final matrix is sorted with respect
% to the concentrations.
%
% Warning: do not use when there are multiple instances of the same time
% point in the time vector (only the first one will be used, and therefore
% an error will be produced).
% 
% Author           : Tjalling Jager 
% Date             : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

n = length(varargin);
t = [];
c = [];
for i = 1:n % run through the supplied matrices
    t_add = varargin{i}(2:end,1);
    if length(unique(t_add)) ~= length(t_add)
        error('The time vector for each data set needs to be unique for mat_combine to work.')
    end
    t = [t;t_add]; % and combine all time vectors
    c = [c varargin{i}(1,2:end)]; % and also all concentration vectors
    N(i) = length(varargin{i}(1,2:end)); % remember length of each concentration vector
end

Ttot = unique(t);

Ctot = nan(length(Ttot)+1,length(c)+1);
Ctot(1,:) = [varargin{1}(1,1) c]; % use the first element of the first dataset for the combined one
Ctot(2:end,1) = Ttot;

locc = 2; % location for columns to put in the data 
for i = 1:n % run through the supplied matrices
    for j = 1:length(Ttot) % run through all time points
        [a,loct] = ismember(Ttot(j), varargin{i}(2:end,1)); % is this timepoint in this matrix?
        if a == 1 % if so, add it in the right location
            Ctot(j+1,locc:locc+N(i)-1) = varargin{i}(loct+1,2:end);
        end
    end
    locc = locc + N(i); % update the location for the next columns 
end

if estim == 1
    % replace NaNs by estimated data ...
    for i = 1:length(c) % run through the concentration vector
        Ctry = Ctot(2:end,i+1); % take the data from the corresponding column
        if any(isnan(Ctry)) % are there NaNs in there?
            loc_nan = 1+find(isnan(Ctry)); % find the NaNs
            loc_ok = 1+find(~isnan(Ctry)); % find the data that are not NaNs
            % and interpolate with cubic spline
            Ctot(loc_nan,i+1) = interp1(Ctot(loc_ok,1),Ctot(loc_ok,i+1),Ctot(loc_nan,1),'pchip','extrap');
        end
    end
end

% sort the columns such that the concentration vector is increasing
Ctot(:,2:end) = sortrows((Ctot(:,2:end))',1)';