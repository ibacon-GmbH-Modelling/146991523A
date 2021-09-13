%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files in the directory 'simple_compound'. As input,
% it gets:
% 
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

% NOTE: glo is now no longer a global here, but passed on in the function
% call from call_deri. That saves 20% calculation time!

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

X = max(X,0); % ensure the states cannot become negative due to numerical issues

Dw = X(glo.locD); % state is the scaled damage (referenced to water)
L  = X(glo.locL); % state is body length
% Rc = X(glo.locR); % state is cumulative reproduction (not used)
S  = X(glo.locS); % state is survival probability

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

% unpack globals
FBV    = glo.FBV;    % dry weight of egg as fraction of dry body weight (-)
KRV    = glo.KRV;    % part. coeff. repro buffer and structure (kg/kg)
kap    = glo.kap;    % approximation for kappa (-)
yP     = glo.yP;     % product of yVA and yAV (-)
Lm_ref = glo.Lm_ref; % reference max length for scaling rate constants
% Note: a reference Lm is needed to properly compare different data sets,
% or when calibrating on more than one data set. If Lm differs, you don't
% want to have different rate constants at the same length.

% unpack model parameters for the basic life history
L0   = par.L0(1);   % body length at start (mm)
Lp   = par.Lp(1);   % body length at puberty (mm)
Lm   = par.Lm(1);   % maximum body length (mm)
rB   = par.rB(1);   % von Bertalanffy growth rate constant (1/d)
Rm   = par.Rm(1);   % maximum reproduction rate (#/d)
f    = par.f(1);    % scaled functional response (-)
hb   = par.hb(1);   % background hazard rate (d-1)

% unpack extra parameters for specific cases
Lf   = par.Lf(1);   % body length at half-saturation feeding (mm)
Lj   = par.Lj(1);   % body length at which acceleration stops (mm)
Tlag = par.Tlag(1); % lag time for start development (d)

% unpack model parameters for the response to toxicants
kd   = par.kd(1);   % dominant rate constant (d-1)
zb   = par.zb(1);   % effect threshold energy budget ([C])
bb   = par.bb(1);   % effect strength energy-budget effects (1/[C])
zs   = par.zs(1);   % effect threshold survival ([C])
bs   = par.bs(1);   % effect strength survival (1/([C] d))

%% Extract correct exposure for THIS time point
% Allow for external concentrations to change over time, either
% continuously, or in steps, or as a static renewal with first-order
% disappearance. For constant exposure, the code in this section is skipped
% (and could also be removed). Note that glo.timevar is used to signal that
% there is a time-varying concentration. This option is set in call_deri.

if glo.timevar(1) == 1 % if we are warned that we have a time-varying concentration ...
    c = read_scen(-1,c,t,glo); % use read_scen to derive actual exposure concentration
    % Note: this has glo as input to the function to save time!
    % the -1 lets read_scen know we are calling from derivatives (so need one c)
end

%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs. This is the
% DEBkiss model, with toxicant effects and starvation module, expressed in
% terms of compound parameters, as presented in the publication (Jager,
% 2020).

L = max(1e-3*L0,L); % make sure that body length is not negative or almost zero (extreme shrinking may do that)
% This should not be needed as shrinking is limited at the bottom of this
% function.

if Lf > 0 % to include feeding limitation for juveniles ...
    f  = f / (1+(Lf^3)/(L^3)); % hyperbolic relationship for f with body volume
    % kd = kd*f; % also reduce dominant rate by same factor? (open for discussion!)
end
if Lj > 0 % to include acceleration until metamorphosis ...
    f = f * min(1,L/Lj); % this implies lower f for L<Lj
end

% calculate stress factor and hazard rate
s = bb*max(0,Dw-zb); % stress level for metabolic effects
h = bs*max(0,Dw-zs); % hazard rate for effects on survival

% Define specific stress factors s*, depending on the mode of action as
% specified in the vector with switches glo.moa.
Si = glo.moa * s;  % vector with specific stress factors from switches for mode of action
sA = min(1,Si(1)); % assimilation/feeding (maximise to 1 to avoid negative values for 1-sA)
sM = Si(2);        % maintenance (somatic and maturity)
sG = Si(3);        % growth costs
sR = Si(4);        % reproduction costs  
sH = Si(5);        % also include hazard to reproduction

% calcululate the actual derivatives, with stress implemented
dL = rB * ((1+sM)/(1+sG)) * (f*Lm*((1-sA)/(1+sM)) - L); % ODE for body length

fR = f; % if there is no starvation, f for reproduction is the standard f
% starvation rules can modify the outputs here
if dL < 0 % then we are looking at starvation and need to correct things
    fR = (f - kap * (L/Lm) * ((1+sM)/(1-sA)))/(1-kap); % new f for reproduction alone
    if fR >= 0  % then we are in the first stage of starvation: 1-kappa branch can help pay maintenance
        dL = 0; % stop growth, but don't shrink
    else        % we are in stage 2 of starvation and need to shrink to pay maintenance
        fR = 0; % nothing left for reproduction
        dL = (rB*(1+sM)/yP) * ((f*Lm/kap)*((1-sA)/(1+sM)) - L); % shrinking rate
    end
end
        
R  = 0; % reproduction rate is zero, unless ... 
if L >= Lp % if we are above the length at puberty, reproduce
    R = max(0,(exp(-sH)*Rm/(1+sR)) * (fR*Lm*(L^2)*(1-sA) - (Lp^3)*(1+sM))/(Lm^3 - Lp^3));
    % Note: hazard to reproduction added with sH
end

dRc = R;             % cumulative reproduction rate
dS  = -(h + hb) * S; % change in survival probability (incl. background mort.)

% For the damage dynamics, there are four feedback factors x* that obtain a
% value based on the settings in the configuration vector glo.feedb: a
% vector with switches for various feedbacks: [surface:volume on uptake,
% surface:volume on elimination, growth dilution, losses with
% reproduction].
Xi = glo.feedb .* [Lm_ref/L,Lm_ref/L,(3/L)*dL,R*FBV*KRV]; % multiply switch factor with feedbacks
xu = Xi(1); % factor for surf:vol scaling uptake rate 
xe = Xi(2); % factor for surf:vol scaling elimination rate 
xG = Xi(3); % factor for growth dilution
xR = Xi(4); % factor for losses with repro

% if switch for surf:vol scaling is zero, the factor must be 1 and not 0!
% Note: this was previously done with a max-to-1 command. However, that is
% not a good idea in combination with Lm_ref (which is not necessarily
% equal to or larger than Lm).
if Xi(1) == 0
    xu = 1;
end
if Xi(2) == 0
    xe = 1;
end

xG = max(0,xG); % stop reverse growth dilution
% NOTE NOTE: reverse growth dilution (concentration by shrinking) is now
% turned OFF as it leads to runaway situations that lead to failure of the
% ODE solvers. However, this needs some further thought!

W   = 1; % TEST optional factor to increase damage (e.g., a relative 'damage partition coeff.')
dDw = kd * (xu * W * c - xe * Dw) - (xG + xR) * Dw; % ODE for scaled damage

if L <= 0.5 * L0 % if an animal has size less than half the start size ...
    dL = 0; % don't let it grow or shrink any further (to avoid numerical issues)
end

dX = zeros(size(X)); % initialise with zeros in right format
if t >= Tlag % when we are past the lag time ...
    dX([glo.locD glo.locL glo.locR glo.locS]) = [dDw;dL;dRc;dS]; % collect all derivatives in one vector dX
end
