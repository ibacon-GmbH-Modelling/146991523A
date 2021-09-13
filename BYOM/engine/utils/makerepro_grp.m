function [data,w] = makerepro_grp(R,S,F,F_fem)

% Usage: [data,w] = makerepro_grp(R,S,F,F_fem)
%
% This is a translation of repro data for grouped animals (data are sum of
% the offspring from all mothers in the replicate), where sex is only
% determined at the end of the test (or not at all).
% 
% Inputs:
%   R: repro matrix with observed offspring at the end of each interval (not cumulated!)
%   S: survivor matrix, with same size as R
%   F: matrix with number of females at the point of sex determination (1 time point per replicate!)
%   F_fem: presumed sex ratio of the animals (used is F is empty and for
%          accounting for mortality)
% 
% Outputs data and w can be included into the cell arrays DATA and W of
% BYOM, at the right place for the corresponding state variable.
%
% Author     : Tjalling Jager 
% Date       : June 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

if ~isequal(size(S),size(R))
    error('The matrices for repro and survivors need to be equal sized')
end

if isempty(F) % then there is no matrix with number of females
    
    % so we'll assume there is a single input for the sex ratio
    F_new = S;
    F_new(2:end,2:end) = F_new(2:end,2:end) * F_fem; % the F_new (number of females) is the survivors times by the sex ratio
    % sex ratio is female:male
    
else % we have final counts of females in the data set
    
    if size(F,2) ~= size(S,2)
        error('The matrix for numbers of females should have the same number of treatments as the repro matrix')
    end
    
    % Create an F_new that is the same size as S and R
    F_new = nan(size(S)); % create new matrix F that is same size as S, filled with NaNs
    F_new(:,1) = S(:,1); % copy first column of S (times)
    F_new(1,2:end) = S(1,2:end); % copy first row of S (replicate IDs)
    [~,ind_F] = ismember(F(2:end,1),F_new(2:end,1)); % locate where the same time points are used
    F_new(ind_F+1,2:end) = F(2:end,2:end); % put the values of F in the right position in F_new
    
    for i = 1:size(F_new,2)-1 % run through replicates
        ind_sex = find(~isnan(F_new(2:end,i+1))); % find where there is an observation
        % At this point, I assume there is only 1 time point at which sex is
        % determined, and that this is at the end of the test (for this
        % replicate)
        if length(ind_sex) ~= 1
            error('Need 1 (and no more than 1) observations on sex')
        end
        F_new(2:ind_sex,i+1) = F_new(ind_sex+1,i+1); % replace previous NaNs with the nr of females
        % if ind_sex+1 < size(F_new,1) % then the sex determination is not at the end of the data set
        if any(isnan(S(2:ind_sex,i+1)))
            error('There can (for now) not be NaNs in the survivor matrix before the sex determination')
        end
        % At this point, I assume there are no NaNs in the survivor matrix
        % before the time at which sex is determined. Otherwise, things become
        % even more complex.
    end
    
    % create a matrix with number of dead females, assuming F_fem of the deaths
    % have been females ...
    F_dead = F_fem*-diff(S(2:end,2:end)); % matrix with estimated dead females (difference in survivors)
    F_dead = [F_dead;zeros(1,size(F_dead,2))]; % append a row with zeros at the end
    F_dead(isnan(F_dead)) = 0; % replace NaNs by zeros
    F_dead = cumsum(F_dead,1,'reverse'); % cumulate back over time ...
    F_new(2:end,2:end) = F_new(2:end,2:end) + F_dead; % add the dead females to the weights matrix
    
end

% Finally, account for the deaths of females within an interval: the
% observed offspring have been produced by the AVERAGE number of females
% alive in the interval (this may also be ignored, then the offspring are
% produced by the females alive AT THE END of the interval). Alternatively,
% we could also use the number alive at the START of the interval ...
F_new(3:end,2:end) = (F_new(3:end,2:end) + F_new(2:end-1,2:end))/2;
w = F_new(2:end,2:end); % the weight matrix does not need the times and replicate IDs

% The cumulated repro per female is calculated by dividing the observed
% offspring production by the estimated average number of females alive,
% and cumulating that over time.
data = R; % start with a copy of R
data(2:end,2:end) = cumsum(R(2:end,2:end) ./ w,1); % new data set as mean cumulative per female!

