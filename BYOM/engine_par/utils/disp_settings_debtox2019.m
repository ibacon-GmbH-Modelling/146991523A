% A small script to display extra information on screen, regarding the
% settings in the DEBtox2019 script. This is useful for archiving, so it is
% always clear what settings were used. Note that this file is also useful
% for extended GUTS analyses (when ODE solver is used). 
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
disp('Extra information on DEBtox2019 (or extended GUTS) fit')
if exist('MOA','var') % if this variable is defined in main script ...
    disp(['Mode of action        : ',num2str(MOA(end,:))])
elseif isfield(glo,'moa')
    disp(['Mode of action        : ',num2str(glo.moa)])
    MOA = [];
else
    MOA = [];
end
if exist('FEEDB','var') % if this variable is not defined in main script ...
    disp(['Feedback configuration: ',num2str(FEEDB(end,:))])
elseif isfield(glo,'feedb')
    disp(['Feedback configuration: ',num2str(glo.feedb)])
    FEEDB = [];    
else
    FEEDB = [];
end

if size(MOA,1)>1 || size(FEEDB,1)>1
    disp('   these are only the settings for the LAST run')
end
disp(' ')
disp('Model and ODE solver ----------------------------------------------')
disp(['   Breaking time vector in call_deri: ',num2str(glo.break_time)])
disp(['   Settings for ODE solver call_deri: ',num2str(glo.stiff)])
if isfield(glo,'Tbp')
    disp(['   Brood-pouch delay for model (d)  : ',num2str(glo.Tbp)])
end
if opt_optim.type == 4 % for the parspace explorer
    disp('Options for parspace explorer -------------------------------------')
    disp(['   Rough settings (recommended: 1)  : ',num2str(opt_optim.ps_rough)])
    disp(['   Make profile likelihoods         : ',num2str(opt_optim.ps_profs)])
    disp('Options for startgrid ---------------------------------------------')
    disp(['   Skip automatic startgrid ranges  : ',num2str(skip_sg)])
elseif opt_optim.type == 1 % for simple optimisation
    disp('Options for simplex optimisation-----------------------------------')
    disp(['   Number of sequential runs        : ',num2str(opt_optim.simno)])
end
disp('===================================================================')

diary off  % close results.out