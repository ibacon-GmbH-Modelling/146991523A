function logprob = calc_prior(par,glo2)

% Usage: logprob = calc_prior(par,glo2)
%
% From structure in <par>, calculate prior probabilities. This
% function does not use the statistics toolbox anymore but calculates the
% probability density for the normal and lognormal distributions directly.
% Only for the beta distribution, the statistics toolbox would be needed.
%
% Author     : Tjalling Jager
% Date       : February 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

namesp = glo2.namesp;
pri    = glo2.pri; % pri will be used in an 'eval' below

logprob = 0;
for i = 1:length(namesp) % run through all parameters that have a prior distribution
    
    if isfield(par,namesp{i}) % is there a parameter matching this prior?
        
        parvals = eval(['par.',namesp{i}]); % extract parameter info to be evaluated for prior prob
        privals = eval(['pri.',namesp{i}]); % extract prior info for calculation of prior prob
        
        % extract parameters for the prior distribution
        a = privals(2);
        b = privals(3);
        
        if privals(1) == 1 % other distributions than beta will work
            if exist('betapdf','file')~=2 % only when betapdf exists as an m-file in the path
                error('You cannot use the beta distribution for priors; you need the statistics toolbox for that.')
            end
        end
        
        switch privals(1) % select correct distribution type
            case 1 % beta distribution
                priprob = betapdf(parvals(1),a,b); % calculate value from PDF
                % a and b are the shape parameters
            case 2 % triangular distribution
                c = privals(4);
                priprob = 0;
                if parvals(1)>= a && parvals(1)<= c
                    priprob = 2*(parvals(1)-a)/((b-a)*(c-a));
                elseif parvals(1)> c && parvals(1)<= b
                    priprob = 2*(b-parvals(1))/((b-a)*(b-c));
                end
                % a is min, b is max, c is center
            case 3 % normal distribution
                % priprob = normpdf(parvals(1),a,b);
                % a is mean, b is standard deviation
                priprob = exp(-0.5 * ((parvals(1) - a)/b)^2) / (sqrt(2*pi)*b);
            case 4 % lognormal distribution
                % priprob = lognpdf(parvals(1),a,b);
                % a is mean of log par, b is standard deviation of log par; base of logarithm is e!
                priprob = exp(-0.5 * ((log(parvals(1)) - a)/b)^2) / (parvals(1) * sqrt(2*pi) * b);
            otherwise
                % do nothing ...
                error('Unknown type of prior distribution.')
        end
        
        priprob = max(1e-100,priprob); % avoid values of zero, as we need to take the log
        logprob = logprob + log(priprob); % add log prob to previous one (sum over all pars)
        
    else
        % error('One of the priors does not have a matching parameter!')
    end
end