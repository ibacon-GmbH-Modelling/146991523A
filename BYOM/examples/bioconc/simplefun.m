%% BYOM function simplefun.m (the model as explicit equations)
%
%  Syntax: Xout = simplefun(t,X0,par,c,glo)
%
% This function calculates the output of the model system. It is linked to
% the script files <byom_bioconc_extra.html byom_bioconc_extra.m> and
% <byom_bioconc_start.html byom_bioconc_start.m>. As input, it gets:
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
% used in this function. The treatments differ in their initial exposure
% concentration that is set in X0mat. Also, the structure _glo_ is not used
% here.

function Xout = simplefun(t,X0,par,c,glo)

%% Unpack initial states
% The state variables enter this function in the vector _X_0.

Cw0 = X0(1); % state 1 is the external concentration at t=0
Ci0 = X0(2); % state 2 is the internal concentration at t=0

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

kd   = par.kd(1);   % degradation rate constant, d-1
ke   = par.ke(1);   % elimination rate constant, d-1
Piw  = par.Piw(1);  % bioconcentration factor, L/kg
Ct   = par.Ct(1);   % threshold external concentration that stops degradation, mg/L

%% Calculate the model output
% This is the actual model, specified as explicit function(s):
%
% $$ C_w=C_{w0} \exp(-k_d t) $$
%
% $$ C_i= a \exp(-k_d t) + (C_{i0}-a) \exp(-k_e t) $$
%
% $$ \textrm{with:} \quad a= \frac{C_{w0} k_e P_{iw}}{k_e-k_d} $$
% 
% *Note:* this function is not completely identical to the solution of the
% system of ODE's. Here, I ignore the stop in degradation, which would
% require quite some fiddling around. Watch out that this solution is only
% valid when _ke_ is NOT equal to _kd_!

Cw = Cw0 * exp(-kd*t); % concentration in water

a  = Cw0 * ke * Piw/(ke-kd); % this factor occurs twice in the solution
Ci = a*exp(-kd*t) + (Ci0 - a) * exp(-ke*t); % internal concentration

Xout = [Cw Ci]; % combine them into a matrix