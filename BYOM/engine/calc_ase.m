function calc_ase(par,opt_ase)

% Usage: calc_ase(par,opt_ase)
% 
% Calculates the asymptotic standard error, coefficient of variation, and
% correlation matrix from the maximum likelihood parameter estimates given
% in total parameter vector <par>.
%
% This function is slightly modified from the function in DEBtoxM.
%
% Basis for the calculation is a numerical derivation of the second
% derivative of the likelihood function to the parameters (the Hessian).
% All derivatives are calculated by taking three points around the ML
% estimate (the estimate itself and +- 1%), and calculate the function
% values. The derivative is then calculated by linear regression on these
% three points.
%
% Note: this is a rather crude estimation of the second derivative, and the
% routine is not very stable: MatLab often has problems to invert the
% matrix. It is MUCH better to use the profile likelihoods instead!

% Author     : Tjalling Jager
% Date       : Nov. 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 DATA
% Note: the calculation requires ttot and ctot to be defined as globals,
% which is done in prelim_checks.m. These parameters are required by
% transfer.m. Only the names parameter is used here for plotting and
% reporting.

% read options from structure
step = opt_ase.step; % relative step size used to calculate derivatives

% extract parameters from the general globals, glo and glo2
n_D   = glo2.n_D;
n_X   = glo2.n_X;
names = glo2.names;

pmat    = packunpack(1,par,0);  % transform structure into a regular matrix
par_sel = pmat(:,2);            % this is the selection vector
pmat(pmat(:,5)==0,1) = log10(pmat(pmat(:,5)==0,1));

% trim parametervectors to use only the reduced set of parameters
p_index  = find(par_sel==1);
parshat  = pmat(p_index,1); % Only the parameters that were fitted

disp(' ')
disp('Calculation of asymptotic standard errors')
disp('Using numerical approximation of the Hessian')
disp('WARNING: this procedure is NOT very accurate, and may even lead to complete nonsense!')
disp('   (especially for small data sets and when there are strong correlations)')
disp('   Profiling the likelihood is MUCH more robust and accurate, also for small data sets')

np   = length(parshat);
Hess = zeros(np,np);

for i = 1:np
    for j = 1:i
        step_i = parshat(i) * step; % take the steps relative to the parameter value
        step_j = parshat(j) * step;
        % Should we catch cases where param estim is (close to) zero??

        parstemp1 = parshat;
        parstemp2 = parshat;
        parstemp3 = parshat;
        parstemp4 = parshat;

        parstemp1(j) = parstemp1(j) + step_j;
        parstemp1(i) = parstemp1(i) + step_i;

        parstemp2(j) = parstemp2(j) - step_j;
        parstemp2(i) = parstemp2(i) + step_i;

        parstemp3(j) = parstemp3(j) + step_j;
        parstemp3(i) = parstemp3(i) - step_i;

        parstemp4(j) = parstemp4(j) - step_j;
        parstemp4(i) = parstemp4(i) - step_i;

        L1 = -1 * transfer(parstemp1,pmat); % calculate LOG likelihood using transfer function
        L2 = -1 * transfer(parstemp2,pmat); % calculate LOG likelihood using transfer function
        L3 = -1 * transfer(parstemp3,pmat); % calculate LOG likelihood using transfer function
        L4 = -1 * transfer(parstemp4,pmat); % calculate LOG likelihood using transfer function

        Hess(i,j) = (L1-L2-L3+L4)/(4*step_i*step_j);
        if i ~= j % all terms are mirrored in the diagonal!
            Hess(j,i) = Hess(i,j);
        end

    end
end

% use them to make covariance matrix, sd's and correlations
Covmat = inv(-1 * Hess);           % covariance matrix
sdvect = sqrt(diag(Covmat));       % standard deviations
CVvect = abs(sdvect ./ parshat);   % coefficient of variation
cormat = Covmat./(sdvect*sdvect'); % correlation matrix

% display the standard errors of the estimates
diary (glo.diary) % collect output in the diary "results.out"
disp(' ')
disp('Estimated parameters, standard errors (CV), and approximate confidence intervals')
disp('=================================================================================')
sdmat = zeros(size(par_sel));
sdmat(p_index) = sdvect;
CVmat = zeros(size(par_sel));
CVmat(p_index) = CVvect;

n = 0;
ndata = n_X * n_D; % total number of data sets
for i = 1:ndata
    n = n + sum(~isnan(DATA{i}(:))); % count number of data points
end
p = 0;
for i = 1:length(names) % run through all parameters
    p = p + eval(['par.',names{i},'(2)']); 
end

if exist('tinv','file') == 2 % then it is a function in the path
    crit = tinv(0.975,n-p); % number of data points minus number of estimated parameters
    disp('Confidence intervals calculated from t-distribution.')
else
    crit = 1.96; % otherwise just use the standard normal distribution
    disp('Confidence intervals calculated from normal distribution (valid when there are a lot of data points).')
end
if any(pmat(:,5)==0)
    disp('Note: some parameters have been estimated on log-scale (shown on log10 scale).')
end

for i = 1:length(p_index)
    Clow = pmat(p_index(i),1) - sdmat(p_index(i))*crit;
    Chi  = pmat(p_index(i),1) + sdmat(p_index(i))*crit;
    fprintf('%-6s %10.4g  se %10.4g (%10.4g) CI %10.4g - %1.4g \n',names{p_index(i)},pmat(p_index(i),1),sdmat(p_index(i)),CVmat(p_index(i)),Clow,Chi)
end

% and display the correlation matrix
disp('=================================================================================')
disp('  ')
disp('Correlation matrix of estimated parameters')
disp('=================================================================================')
fprintf('%7s','  ')
for i = 1:length(p_index)
    fprintf('%11s ',names{p_index(i)})
end
fprintf('\n')
for i = 1:length(p_index)
   
    fprintf('%-6s ',names{p_index(i)})
    for j = 1:i % length(p_index),
        fprintf('%11.4g ',cormat(i,j))
    end
    fprintf('\n')
end
if any(pmat(:,5)==0)
    disp('Note: some parameters have been estimated on log-scale (shown on log10 scale).')
end
disp('=================================================================================')
disp(['Time required: ' secs2hms(toc)])
diary off