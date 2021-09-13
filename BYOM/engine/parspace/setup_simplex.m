function [phat,fval,exitflag] = setup_simplex(pfit,option,varargin)

% Usage: [phat,fval,exitflag] = setup_simplex(pfit,option,varargin)
% 
% A wrapper function around the optimisation routine. This allows to change
% the optimisation routine without needing to change all the calls in the
% rest of the software. However, in the BYOM version, only Matlab's
% <fminsearch> is implemented. The main reason to keep this wrapper is to
% stay close to the openGUTS version of the code in the other parspace
% functions.
% 
% Inputs
% <pfit>        vector with parameters to be fitted (initial values)
% <option>      1) standard optimisation, 0) rough optimisation
% <varargin>    additional arguments to be passed on to <transfer>
% 
% Outputs
% <phat>        vector with optimised parameter values
% <fval>        optimised value of minus-log-likelihood
% <exitflag>    1) convergence, 0) no convergence or errors
% 
% Author     : Tjalling Jager
% Date       : May 2020
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl).  This file
% is a slightly modified version of the <setup_simplex> code that is
% distributed as part of the Matlab version of openGUTS (see
% http://www.openguts.info). Therefore, this code is distributed under the
% same license as openGUTS (GPLv3). The modifications are only to ensure
% that the code operates in the general BYOM framework.
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%  
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

% Always use <fminsearch> in BYOM, no need for <nelmin>.

% Make options for the standard and rough simplex optimisations.
oldopts = optimset('fminsearch'); % retrieve default options for <fminsearch>
opts    = optimset(oldopts,'FunValCheck','on','Display','off'); % function check on, show no output at all
if option == 0 % options for rough simplex
    opts = optimset(opts,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',30*length(pfit)); % rough Simplex options, based on number of fitted parameters
    % Note: these rough options will use a smaller total number of function
    % evaluations: quick and dirty.
end

% Use <fminsearch> and the standard <transfer> function to fit.
[phat,fval,exitflag] = fminsearch('transfer',pfit,opts,varargin{:}); % do a normal optimisation

