function [par_out,error_flag] = plot_grid_ind(filename)

% Usage: [par_out,error_flag] = plot_grid_ind(filename)
% 
% Function to combine MAT files for two chemicals (as generated by the
% parameter-space explorer) to make a MAT file for the mixture, assuming
% independent action. Independent action (in the strictest sense) implies
% that the two chemicals do not interact in any way. Therefore, combining
% multiple parameter clouds is simply a matter of randomly combining sets
% from all clouds. The only important thing is to make sure that the
% resulting cloud is not too large. This function only works properly if
% the MAT file was generated by the BYOM GUTS package for the standard
% reduced models! Output is also meant for the standard GUTS binary mixture
% package (which will be released at some point).
% 
% filename   cell array with strings for the name of MAT files to combine
% 
% Author     : Tjalling Jager
% Date       : July 2021
% Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo2

SETTINGS_OPTIM = setup_settings(1); % load settings from file (use rough settings)
names = {'kd';'mw';'hb';'bw';'Fs'}; % names of the model parameters in the saved mat files
% this assumes that the standard GUTS package was used to create these files!

n_s = length(filename); % number of samples to compare
if n_s > 4 % for now, we only allow two files anyway, so this is for the future
    error_str = ('Cannot compare more than 4 samples at the moment');
    f = errordlg(error_str,'Error no. of sets'); % display a message box with the errors
    error_flag = 1; % signal main script that we ran into errors
    return % simply stop as there is no useful input to work with
end

%% Run through the samples to load files and collect some things

for i_s = 1:n_s % run through samples
    
    % Load sample <i_s>
    if exist([filename{i_s}],'file') == 2 % check if it exists first.
        load([filename{i_s}],'pmat','coll_all','pmat_print')
        % Load all of the saved information from the <calc_parspace> run.
    else % otherwise, produce an error (I don't think that is possible anymore)
        error(['There is no confidence set with filename ',filename{i_s},' in your working directory, so run calibration (with correct settings) first.'])
    end
    
    sel(i_s) = str2num(filename{i_s}(end-7)); % SD or IT (is part of the file name)
    
    % Extract useful information from <pmat>.
    ind_fit  = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
    n_fit    = length(ind_fit); % number of fitted parameters
    pmat_lim = pmat(ind_fit,:); % limit <pmat_lim> to fitted parameters
    ind_log  = find(pmat_lim(:,5) == 0); % indices to log-scale parameters in fitted <pmat> (vector)
    
    pmat_lim(ind_log,1) = 10.^pmat_lim(ind_log,1); % put back on normal scale, where needed
    coll_all(:,ind_log) = 10.^coll_all(:,ind_log); % also in the sample
    names_lim           = names(ind_fit); % only keep names for fitted parameters
        
    % collect the results for further processing
    BNDS{i_s}  = pmat(ind_fit,[3 4]); % collect bounds used for original optimisation
    COLL{i_s}  = coll_all;
    NAMES{i_s} = names_lim;
    
    ind_hb  = find(strcmp(names,'hb')==1); % find location of <hb> in <pmat>
    HB{i_s} = pmat(ind_hb,1); % collect value for <hb>
    if pmat(ind_hb,2) == 1
        warning('off','backtrace')
        warning(['Background mortality was fitted in dataset ',filename{i_s},', this is ignored!'])
        warning('on','backtrace')
    end
    if sel(i_s)~=sel(1)
        error('all mat files for comparison need to be made with the same death mechanism')
    end
    
end

% Chi2-criteria and indices to parameter sets to plot for inner and outer rim.
% We work here with the log-likelihood itself, so the chi2 criterion needs to be divided by 2.
chicrit_joint  = 0.5 * SETTINGS_OPTIM.crit_table(6,1); % criterion for joint 95% CI or parameters
chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1); % criterion for single-parameter CIs
chicrit_prop   = 0.5 * SETTINGS_OPTIM.crit_prop(2);    % criterion for upper band of propagation band
chicrit_prop2  = 0.5 * SETTINGS_OPTIM.crit_prop(1);    % criterion for lower band of propagation band

%% Create a new <coll_all> matrix for use with mixtures

coll_all = []; % start with empty <coll_all> matrix

if sel(1) == 1 % for SD
    
    % find indices for all model parameters in names vector for each chemical
    ind_kdA = find(strcmp(NAMES{1},'kd')==1);
    ind_kdB = find(strcmp(NAMES{2},'kd')==1);
    ind_mwA = find(strcmp(NAMES{1},'mw')==1);
    ind_bwA = find(strcmp(NAMES{1},'bw')==1);
    ind_mwB = find(strcmp(NAMES{2},'mw')==1);
    ind_bwB = find(strcmp(NAMES{2},'bw')==1);
    
    MLL = COLL{1}(1,end)+COLL{2}(1,end); % best minloglik is simply the sum of the first elements of the final column
    
    % make sure best combination is in the <coll_add> matrix
    coll_add = nan(1,7); % what we're adding to coll_all for set iA
    coll_add(:,1) = COLL{1}(1,ind_kdA); % collect kdA
    coll_add(:,2) = COLL{1}(1,ind_mwA); % collect mwA
    coll_add(:,3) = COLL{1}(1,ind_bwA); % collect bwA
    coll_add(:,4) = COLL{2}(1,ind_kdB); % collect kdB
    coll_add(:,5) = COLL{2}(1,ind_mwA); % collect mwB
    coll_add(:,6) = COLL{2}(1,ind_bwA); % collect bwB
    coll_add(:,7) = COLL{1}(1,end)+COLL{2}(1,end); % new MLL
    
    coll_all = cat(1,coll_all,coll_add); % collect in total <coll_all>
    min_sz   = min(size(COLL{1},1),size(COLL{2},1)); % minimum size of the two clouds
    ind_prop = 0;
    
    f = waitbar(0,'Combining MAT files for SD. Please wait.','Name','plot_grid_ind.m');
    while ind_prop < 5000 % just a few rounds, until we have sufficient sets in the inner rim (incl. propagation set)
        % I aim here for at least 5000 elements in the inner rim of the new
        % cloud, which seems large enough.
        waitbar(max(1,ind_prop/5000),f); % update waiting bar
        
        % random permutations of parameter sets for A and B. <min_sz> elements
        indA = randperm(min_sz);
        indB = randperm(min_sz);
        
        coll_add = nan(min_sz,7); % what we're adding to <coll_all> in this round
        coll_add(:,1) = COLL{1}(indA,ind_kdA); % collect kdA
        coll_add(:,2) = COLL{1}(indA,ind_mwA); % collect mwA
        coll_add(:,3) = COLL{1}(indA,ind_bwA); % collect bwA
        coll_add(:,4) = COLL{2}(indB,ind_kdB); % collect kdB
        coll_add(:,5) = COLL{2}(indB,ind_mwB); % collect mwB
        coll_add(:,6) = COLL{2}(indB,ind_bwB); % collect bwB
        coll_add(:,7) = COLL{1}(indA,end)+COLL{2}(indB,end); % minloglik is sum of last columns
        
        coll_add(coll_add(:,7)>MLL+chicrit_joint,:) = []; % remove ones that are pretty bad
        
        coll_add  = sortrows(coll_add,size(coll_add,2)); % sort what we're adding based on minloglik
        ind_fin95 = find(coll_add(:,end) < coll_add(1,end) + chicrit_joint,1,'last'); % index to last element of sample that is still in joint CI
        ind_prop  = find(coll_add(:,end) < coll_add(1,end) + chicrit_prop,1,'last');  % index to last element of sample that is still in propagation band

        % we don't need so many in the outer rim, so we can downsample
        thin_sz  = ind_fin95-ind_prop; % elements in outer rim (outside propagation set)
        ind_thin = randperm(thin_sz,min(thin_sz,ind_prop)); % downsample to minimum of inner rim and outer rim
        % This means that you don't take more points from the outer rim
        % than there are in the inner. That seems to be okay.
        
        coll_all = cat(1,coll_all,coll_add(1:ind_prop,:),coll_add(ind_prop+ind_thin,:)); % collect in total <coll_all>
        coll_all = sortrows(coll_all,size(coll_all,2)); % sort again on minloglik
        ind_prop = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop,1,'last'); % index to last element of sample that is still in propagation band
        
    end
    close(f) % close the waiting bar
    
elseif sel(1) == 2 % for IT

    % find indices for all model parameters in names vector for each chemical
    ind_kdA = find(strcmp(NAMES{1},'kd')==1);
    ind_kdB = find(strcmp(NAMES{2},'kd')==1);
    ind_mwA = find(strcmp(NAMES{1},'mw')==1);
    ind_FsA = find(strcmp(NAMES{1},'Fs')==1);
    ind_mwB = find(strcmp(NAMES{2},'mw')==1);
    ind_FsB = find(strcmp(NAMES{2},'Fs')==1);
    
    MLL = COLL{1}(1,end)+COLL{2}(1,end); % best minloglik is simply the sum of the first elements of the final column
    
    % make sure best combination is in the <coll_add> matrix
    coll_add = nan(1,7); % what we're adding to <coll_all> in this round
    coll_add(:,1) = COLL{1}(1,ind_kdA); % collect kdA
    coll_add(:,2) = COLL{1}(1,ind_mwA); % collect mwA
    coll_add(:,3) = COLL{1}(1,ind_FsA); % collect FsA
    coll_add(:,4) = COLL{2}(1,ind_kdB); % collect kdB
    coll_add(:,5) = COLL{2}(1,ind_mwA); % collect mwB
    coll_add(:,6) = COLL{2}(1,ind_FsA); % collect FsB
    coll_add(:,7) = COLL{1}(1,end)+COLL{2}(1,end); % minloglik is sum of last columns
    
    coll_all = cat(1,coll_all,coll_add); % collect in total <coll_all>
    
    min_sz   = min(size(COLL{1},1),size(COLL{2},1)); % minimum size of the two clouds
    ind_prop = 0;
    
    f = waitbar(0,'Combining MAT files for IT. Please wait.','Name','plot_grid_ind.m');
    while ind_prop < 5000 % just a few rounds, until we have sufficient sets in the inner rim (incl. propagation set)
        % I aim here for at least 5000 elements in the propagation bound of
        % the new cloud, which seems large enough.
        waitbar(max(1,ind_prop/5000),f); % update waiting bar
        
        % random permutations of parameter sets for A and B. <min_sz> elements
        indA = randperm(min_sz);
        indB = randperm(min_sz);
        
        coll_add = nan(min_sz,7); % what we're adding to <coll_all> in this round
        coll_add(:,1) = COLL{1}(indA,ind_kdA); % collect kdA
        coll_add(:,2) = COLL{1}(indA,ind_mwA); % collect mwA
        coll_add(:,3) = COLL{1}(indA,ind_FsA); % collect FsA
        coll_add(:,4) = COLL{2}(indB,ind_kdB); % collect kdB
        coll_add(:,5) = COLL{2}(indB,ind_mwB); % collect mwB
        coll_add(:,6) = COLL{2}(indB,ind_FsB); % collect FsB
        coll_add(:,7) = COLL{1}(indA,end)+COLL{2}(indB,end); % minloglik is sum of last columns
        
        coll_add(coll_add(:,7)>MLL+chicrit_joint,:) = []; % remove ones that are pretty bad
        
        coll_add  = sortrows(coll_add,size(coll_add,2)); % sort what we're adding based on minloglik
        ind_fin95 = find(coll_add(:,end) < coll_add(1,end) + chicrit_joint,1,'last'); % index to last element of sample that is still in joint CI
        ind_prop  = find(coll_add(:,end) < coll_add(1,end) + chicrit_prop,1,'last');  % index to last element of sample that is still in propagation band

        % we don't need so many in the outer rim, so we can downsample
        thin_sz  = ind_fin95-ind_prop; % elements in outer rim (outside propagation set)
        ind_thin = randperm(thin_sz,min(thin_sz,ind_prop)); % downsample to minimum of inner rim and outer rim
        % This means that you don't take more points from the outer rim
        % than there are in the inner. That seems to be okay.
        
        coll_all = cat(1,coll_all,coll_add(1:ind_prop,:),coll_add(ind_prop+ind_thin,:)); % collect in total <coll_all>
        coll_all = sortrows(coll_all,size(coll_all,2)); % sort again based on minloglik
        ind_prop = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop,1,'last'); % index to last element of sample that is still in propagation band
    end
    close(f) % close the waiting bar
    
end

coll_all = sortrows(coll_all,size(coll_all,2)); % sort again based on minloglik
% Find indices to sets within inner and outer rim.
ind_single = find(coll_all(:,end) < coll_all(1,end) + chicrit_single,1,'last'); % index to last element of sample that is still in inner rim
ind_fin95  = find(coll_all(:,end) < coll_all(1,end) + chicrit_joint,1,'last');  % index to last element of sample that is still in joint CI
ind_prop   = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop,1,'last');   % index to last element of sample that is still in upper edge propagation band
ind_prop2  = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop2,1,'last');  % index to last element of sample that is still in lower edge propagation band

%% Create new <pmat> for the mixture

glo2.names = {'hb';'kdA';'mwA';'bwA';'FsA';'kdB';'mwB';'bwB';'FsB';'IAB'}; % names of the model parameters for new coll_all and pmat

pmat      = zeros(10,5); % create a new <pma>t for the mixture independent case
pmat(1,:) = [mean([HB{1} HB{2}])  0    0    1 1]; % for now, take <hb> as mean for two data sets!
if pmat(1,1) > 1e-6 % best to warn that this is happening
    warning('off','backtrace')
    warning('Background hazard hb is set to MEAN for two data sets!')
    warning('on','backtrace')
end

% fill <pmat> with information from new <coll_all> for parameters common to SD and IT
pmat(2,:)  = [coll_all(1,1) 1 min(coll_all(:,1)) max(coll_all(:,1)) 0]; % kdA
pmat(3,:)  = [coll_all(1,2) 1 min(coll_all(:,2)) max(coll_all(:,2)) 0]; % mwA
pmat(6,:)  = [coll_all(1,4) 1 min(coll_all(:,4)) max(coll_all(:,4)) 0]; % kdB
pmat(7,:)  = [coll_all(1,5) 1 min(coll_all(:,5)) max(coll_all(:,5)) 0]; % mwB
pmat(10,:) = [0     0 -100  100 1]; % IAB

% Make the bounds a bit wider, but not wider than the original bounds in the saved pmat
bnds_extra = 1.5; % what to divide/multiply min-max bounds with to make them a bit wider than sample
% set same bounds for kdA and kdB, and for mwA and mwB
pmat([2 6],3) = max(pmat([2 6],3)/bnds_extra,[BNDS{1}(ind_kdA,1);BNDS{2}(ind_kdB,1)]); % min bounds for kdA and kdB
pmat([2 6],4) = min(pmat([2 6],4)*bnds_extra,[BNDS{1}(ind_kdA,2);BNDS{2}(ind_kdB,2)]); % max bounds for kdA and kdB
pmat([3 7],3) = max(pmat([3 7],3)/bnds_extra,[BNDS{1}(ind_mwA,1);BNDS{2}(ind_mwB,1)]); % min bounds for mwA and mwB
pmat([3 7],4) = min(pmat([3 7],4)*bnds_extra,[BNDS{1}(ind_mwA,2);BNDS{2}(ind_mwB,2)]); % max bounds for mwA and mwB

% fill <pmat> with information from new <coll_all> for parameters NOT common to SD and IT
if sel(1) == 1 % for SD
    pmat(4,:) = [coll_all(1,3) 1 min(coll_all(:,3)) max(coll_all(:,3)) 0]; % bwA
    pmat(8,:) = [coll_all(1,6) 1 min(coll_all(:,6)) max(coll_all(:,6)) 0]; % bwB
    pmat(5,:) = [3     0    1  20 1]; % FsA
    pmat(9,:) = [3     0    1  20 1]; % FsB
    
    pmat([4 8],3) = max(pmat([4 8],3)/bnds_extra,[BNDS{1}(ind_bwA,1);BNDS{2}(ind_bwB,1)]); % min bounds for bwA and bwB
    pmat([4 8],4) = min(pmat([4 8],4)*bnds_extra,[BNDS{1}(ind_bwA,2);BNDS{2}(ind_bwB,2)]); % max bounds for bwA and bwB
    
elseif sel(1) == 2 % for IT
    pmat(5,:) = [coll_all(1,3) 1 min(coll_all(:,3)) max(coll_all(:,3)) 0]; % FsA
    pmat(9,:) = [coll_all(1,6) 1 min(coll_all(:,6)) max(coll_all(:,6)) 0]; % FsB
    pmat(4,:) = [1e3 0 1e-6 1e6 1]; % bwA
    pmat(8,:) = [1e3 0 1e-6 1e6 1]; % bwB
    
    pmat([5 9],3) = max(pmat([5 9],3)/bnds_extra,[BNDS{1}(ind_FsA,1);BNDS{2}(ind_FsB,1)]); % min bounds for FsA and FsB
    pmat([5 9],4) = min(pmat([5 9],4)*bnds_extra,[BNDS{1}(ind_FsA,2);BNDS{2}(ind_FsB,2)]); % max bounds for FsA and FsB

end

ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)

% bounds that are smaller than a factor 10 can be fitted on normal scale
chk_bnds = pmat(:,2) == 1 & (pmat(:,4)./pmat(:,3)) < 10; % which elements are smaller than 10?
pmat(chk_bnds,5) = 1; % put them to normal scale

%% Display, plot and final things

WRAP.glo  = [];
WRAP.glo2 = glo2;

% display par in a formatted way so they can be directly copied into the code of the main script
par_out = packunpack(2,[],pmat,WRAP);
disp(' ')
if sel(1) == 1
    disp('  Independent action for stochastic death (SD)')
else
    disp('  Independent action for individual tolerance (IT)')
end
disp('You can copy-paste following lines into script for mixture predictions')
disp('if you like to run the parameter-space explorer on the mixture data.')
print_par(par_out) % this prints out the optimised parameter values in a

pmat_lim = pmat(ind_fit,:);
ind_log  = find(pmat_lim(:,5) == 0); % indices to log-scale parameters in fitted <pmat> (vector)

% improvise a <pmat_print>, though without any CIs, so it can be saved
pmat_print = zeros(length(ind_fit),8);
pmat_print(:,1) = pmat_lim(:,1); % copy best-fit values
% add CIs for the model parameters, estimated as the edges of the
% propagation set; this is a good approximation
pmat_print(:,2) = (min(coll_all(1:ind_prop,1:end-1),[],1))';
pmat_print(:,3) = (max(coll_all(1:ind_prop,1:end-1),[],1))';

disp(' ')
disp('=================================================================================')
disp('Results of the parameter-space combination for mixture predictions')
disp('=================================================================================')
if sel(1) == 1
    disp('Independent action for stochastic death (SD)')
else
    disp('Independent action for individual tolerance (IT)')
end
disp(['   Chemical A is taken from: ',filename{1}])
disp(['   Chemical B is taken from: ',filename{2}])
disp(['   Sample: ',num2str(ind_fin95),' sets in joint CI and ',num2str(ind_single),' in inner CI.'])
disp(['   Propagation set: ',num2str(1+ind_prop - ind_prop2),' sets will be used for error propagation.'])

FVAL = coll_all(1,end);
fprintf('   Minus log-likelihood has reached the value %#1.2f (AIC=%#1.2f). \n',FVAL,2*length(ind_fit)+2*FVAL)

disp(['Approx. best estimates and 95% CIs on fitted parameters'])
disp_pmat(pmat,pmat_print,WRAP); % call dedicated function for printing the <pmat>
disp('Confidence intervals are estimated from the sample (inner rim, incl. propagation band)')
disp('so they will generally slightly exaggerate the true intervals.')

% put <pmat> on log-scale for log-parameters (this is how BYOM uses it)
pmat(pmat(:,2)==1 & pmat(:,5)==0,1) = log10(pmat(pmat(:,2)==1 & pmat(:,5)==0,1));

% put <coll_all> on log-scale for log-parameters
coll_all(:,ind_log) = log10(coll_all(:,ind_log));

for i_p = 1:length(ind_fit) % trick to make sure profiles are plotted (but without refinement red line)
    coll_prof_pruned{i_p} = [NaN NaN];
end

% And make a plot of the final results as parameter-space plot
figh = plot_grid(pmat,[],coll_all,coll_prof_pruned,[],SETTINGS_OPTIM,WRAP); 

% user can give a filename for saving, but I try to provide a good default
filenm = ['IND_',filename{1}(11:end-12),'+',filename{2}(11:end)];

% % Instead of the line above, feel free to uncomment the code below to use
% % the Matlab GUI for entering a filename for saving.
% filenm = uiputfile(['IND_',filename{1}(11:end-12),'+',filename{2}(11:end)],'Save new project file for mixture of chemical A and B'); % use Matlab GUI to select name for saving MAT file
% if numel(filenm) == 1 && filenm(1) == 0 % if cancel is pressed ...
%     error_flag = 1; % signal main script that we ran into errors
%     return % simply stop as the user does not want to save
% end

saveplt = 0; % for now, don't save (since this function is not normally called from a standard byom script, the glo.saveplot will not be defined)
figure(figh)  % use existing handle and make sure that the multiplot is the current plot
% Add a title to the plot.
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
text(0.5, 1,['File: ',filenm,' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
if saveplt > 0
    save_plot(figh,['parspace_indmixpred_',glo.basenm]) % save parspace plot in output folder
end
snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output

% save new mat file
save(filenm,'pmat','coll_all','pmat_print','coll_prof_pruned')
