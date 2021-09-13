% A small script to check the data sets entered in the immobility script.
% Only the most obvious errors are spotted!
%
% Author     : Tjalling Jager 
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

data_coll = zeros(size(DATA{3}(3:end,2:end))); % intialise matrix with zeros
for i = 3:6 % run through data sets
    data = DATA{i}(3:end,2:end); % take out relavant data (not first row, not second row, not first column)
    data(isnan(data)) = 0; % replace NaNs by zeros
    
    if ~isequal(DATA{i}(2,2:end),DATA{3}(2,2:end))
        error('The data format is not correct: the second rows in each data set must match (start nr. of individuals in each replicate)')
    end
    if ~isequal(size(DATA{i}),size(DATA{3}))
        error('The data format is not correct: all matrices must be the same size')
    end
    
    data_coll = data_coll + data; % add the data set to the previous one
    
end

for j =1:size(data_coll,1) % run through rows (time points)
    if ~isequal(data_coll(j,:),DATA{3}(2,2:end))
        error('The data format is not correct: the sum across all categories must match the starting number of animals')
    end
end