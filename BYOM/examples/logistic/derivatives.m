%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files byom_logistic_*. The model is for logistic
% population growth, including an optional constant predation rate. As
% input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
% 
% Note: _glo_ is now an input to this function rather than a global. This
% makes the code considerably faster.
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start
% Note that, in this example, variable _c_ (the scenario identifier) is not
% used in this function. The treatments differ in their initial population
% size that is set in X0mat. Also, the time _t_ and structure _glo_ are not
% used here.

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

N = X(1);     % state 1: population size
N = max(N,0); % make sure population size does not become negative

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file. The
% 1 between parentheses is needed for fitting the model, as each parameter
% has 4-5 associated values.

r   = par.r(1);   % growth rate
K   = par.K(1);   % carrying capacity
p   = par.p(1);   % predation rate

%% Calculate the derivatives
% This is the actual model, specified as ODE.

dN = r*N*(1-N/K)-p; % logistic growth with predation

if N == 0
    dN = 0; % a population without individuals should not grow or shrink
end

dX = [dN]; % collect derivative in the vector dX