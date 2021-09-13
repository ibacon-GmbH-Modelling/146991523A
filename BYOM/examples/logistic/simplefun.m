%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the model system. It is linked to
% the script files byom_logistic_*. The model is for logistic population
% growth. Note that the constant predation rate of the ODE version is NOT
% used here. As input, it gets:
%
% * _t_   is the time vector
% * _X0_  is a vector with the initial values for states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
% 
% Note: _glo_ is now an input to this function rather than a global. This
% makes the code faster (though only very little when using simplefun).
%
% Time _t_ is handed over as a vector, and scenario name _c_ as single
% number, by <call_deri.html call_deri.m> (you do not have to use them in
% this function). Output _Xout_ (as matrix) provides the output for each
% state at each _t_.
%
% * Author: Tjalling Jager 
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start
% Note that, in this example, variable _c_ (the scenario identifier) is not
% used in this function. The treatments differ in their initial population
% size that is set in X0mat. Also, the structure _glo_ is not used here.

function Xout = simplefun(t,X0,par,c,glo)

%% Unpack initial states
% The state variables enter this function in the vector _X_0.

N0 = X0(1); % state 1 is the population size at t=0

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

r   = par.r(1);   % growth rate
K   = par.K(1);   % carrying capacity
% p   = par.p(1);   % predation rate (not used in this analytical solution)

%% Calculate the model output
% This is the actual model, specified as explicit function(s). Note that
% this is the model *without* predation. This has a nice solution, while I
% am unsure whether there is a solution for the case with predation (I
% guess that needs Maple).

N = (K * N0 * exp(r*t)) ./(K + N0 * (exp(r*t)-1));

Xout = [N]; % combine them into a matrix