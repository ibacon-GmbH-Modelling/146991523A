function logprob = calc_zerovar(zvd,glo2)

% Usage: logprob = calc_zerovar(zvd,glo2)
%
% From structure in zvd, calculate probabilities of zero-variate data. This
% function does not use the statistics toolbox anymore but calculates the
% probability density for the normal distribution directly.
% 
% Note: zvd is used in an 'eval' setting below.
%
% Author     : Tjalling Jager
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

namesz = glo2.namesz;

logprob = 0;
for i = 1:length(namesz) % run through all parameters that have a zero-variate data point
    
        zvd_mean = eval(['zvd.',namesz{i},'(1)']); % extract mean for this data point (as defined in script)
        zvd_sd   = eval(['zvd.',namesz{i},'(2)']); % extract sd for this data point (as defined in script)
        zvd_eval = eval(['zvd.',namesz{i},'(3)']); % extract calculated value for this data point (as done in call_deri)
        
        % calculate normal probability for this data point
        % zvdprob = normpdf(zvd_eval,zvd_mean,zvd_sd); % this requires the statistics toolbox
        zvdprob = exp(-0.5 * ((zvd_eval - zvd_mean)/zvd_sd)^2) / (sqrt(2*pi)*zvd_sd);
        
        zvdprob = max(1e-200,zvdprob); % avoid values of zero, as we need to take the log
        logprob = logprob + log(zvdprob); % add log prob to previous one (sum over all pars)
        
end