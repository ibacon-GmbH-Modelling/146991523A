% A small script to display extra information on screen, regarding the
% settings in the immobility script. This is useful for archiving, so it is
% always clear what settings were used.
%
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

diary (glo.diary) % collect screen output in the diary "results.out"

disp(' ')
disp('===================================================================')
disp('Extra information on GUTS immobility fit')

if isfield(glo,'sel')
    if glo.sel == 1
        disp('   Death mechanism       : SD')
    elseif glo.sel == 2
        disp('   Death mechanism       : IT')
    end
else
    if exist('SEL','var') % if this variable is defined in main script ...
        if SEL(end,1) == 1
            disp('   Death mechanism (last plot) : SD')
        elseif SEL(end,1) == 2
            disp('   Death mechanism (last plot) : IT')
        end
    end
end
disp('Model and ODE solver ----------------------------------------------')

if isfield(glo,'damconfig')
    disp(['   Damage configuration             : ',num2str(glo.damconfig)])
elseif exist('CONFIG','var') % if this variable is defined in main script ...
    disp(['   Damage configuration (last plot) : ',num2str(CONFIG(end,:))])
end

disp(['   Fast kinetics settings           : ',num2str(glo.fastslow)])
disp(['   Breaking time vector in call_deri: ',num2str(glo.break_time)])
disp(['   Settings for ODE solver call_deri: ',num2str(glo.stiff)])
if opt_optim.type == 4 % for the parspace explorer
    disp('Options for parspace explorer -------------------------------------')
    disp(['   Rough settings (recommended: 1)  : ',num2str(opt_optim.ps_rough)])
    disp(['   Make profile likelihoods         : ',num2str(opt_optim.ps_profs)])
    disp('Options for startgrid ---------------------------------------------')
    disp(['   Skip automatic startgrid ranges  : ',num2str(skip_sg)])
    disp(['   Limit rate constants ke and kr   : ',num2str(glo.lim_k)])
    disp(['   Fit all thresholds on log scale  : ',num2str(glo.mw_log)])
elseif opt_optim.type == 1 % for simple optimisation
    disp('Options for simplex optimisation-----------------------------------')
    disp(['   Number of sequential runs        : ',num2str(opt_optim.simno)])
end
disp('===================================================================')

diary off  % close results.out