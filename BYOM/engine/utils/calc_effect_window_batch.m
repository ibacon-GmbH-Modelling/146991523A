function [file_out,twin] = calc_effect_window_batch(par_plot,Twin,opt_ecx,opt_conf)

% Usage: [file_out,twin] = calc_effect_window_batch(par_plot,Twin,opt_ecx,opt_conf)
%           
% Batch calculations for toxic effects due to an exposure profile, with
% moving time window and different (fixed) multiplication factors. This
% functions uses the Matlab GUI element to select files, and calls
% calc_effect_window repeatedly for the actual calculations. 
% 
% <par_plot>   parameter structure for the best-fit curve; if left empty the
%              structure from the saved sample is used
% <Twin>       length of time window (days)
% <opt_ecx>    options structure for ECx and EPx calculations
% <opt_conf>   options structure for making confidence intervals (needed if
%              par_plot is to be read from file)
% 
% Output <file_out> is a cell array with the filenames of the profiles that
% require a closer look. If the profiles are not located in the working
% directory, <file_out> will contain the entire path to these files. Output
% <twin> will contain the start time for the window that was flagged in
% <file_out>. This is NOT necessarily the window that will yield the lowest
% EPx, though.
% 
% Author     : Tjalling Jager 
% Date       : August 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo

% par_read  = opt_ecx.par_read;  % when set to 1 read parameters from saved set, but do NOT make CIs
MF        = opt_ecx.mf_range;  % range for MFs to make plots
mf_crit   = opt_ecx.mf_crit;   % MF trigger for flagging potentially critical profiles

%% Select and load profiles from text
% These may be located in a different location than the current working
% directory. Therefore, filepath is also used.

if isfield(glo,'moa') && isfield(glo,'feedb')
    if glo.feedb(2) == 1 && any(glo.moa(1:3)==1)
        warning('off','backtrace')
        warning('Calculating effect windows at fixed MFs is NOT recommended for this combination of pMoA and feedbacks!')
        warning('There is a possibility of (limited) higher effects at lower MFs, so use extra care.')
        disp(' '), warning('on','backtrace')
    end
end

% Use Matlab GUI element to select files for loading
[filename,filepath] = uigetfile('input_data/*.txt','Select file(s) with exposure profiles for batch analysis','MultiSelect','on'); % use Matlab GUI to select file(s)
if ~iscell(filename) && numel(filename) == 1 && filename(1) == 0 % if cancel is pressed ...
    return % simply stop
end
if ~iscell(filename) % the GUI makes it a cell array only when multiple files are selected ...
    filename = {filename}; % but we also want a cell array if it's only one
end

% load the best parameter set from file
if isempty(par_plot) % also if par_plot is not provided
    % Then we want to read parameters from file
    if isempty(opt_conf)
        error('opt_conf needs to be provided to read parameter structure from MAT file')
    end
    [~,par]  = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    par_plot = par; % simply use the one from the sample file
    disp('The parameter structure is read from the MAT file')
    opt_ecx.par_read = 0;  % when set to 1 read parameters from saved set, but do NOT make CIs
    % this makes sure that calc_effect_window does not re-read the
    % parameters every time
end
if ~isempty(opt_conf)
    opt_conf.type = 0; % use values from slice sampler (1), likelihood region(2) to make intervals
    % for batch calculations we do not make CIs!
end

%% Run through profiles
% This calls calc_effect_window every time, for each profile.

opt_ecx.batch_epx = 1; % when set to 1 use batch mode (no output to screen)

COLL = [];
% If the profile is not located in the working directory, return the partial path
if ~strcmp(filepath(1:end-1),pwd)
    a = strfind(filepath,filesep);
    if length(a) > 2
        disp(['profiles from directory: ',filepath(a(end-2):end)])
    else
        disp(['profiles from directory: ',filepath])
    end
end

f = waitbar(0,'Calculating moving time window in batch mode.','Name','calc_effect_window_batch.m');

for i = 1:length(filename)
    
    disp(['Calculating profile from file: ',filename{i},' (',num2str(i),' of ',num2str(length(filename)),')'])
        
    % <MinColl> collects, for each state, the MF, minimum of the trait relative
    % to the control (thus largest effect), and time at which it occurs.
    % <MinCI> collects the CI for the minimum trait level. Structure is
    % MinColl{i_MF}(i_X,:). Output of ind_traits is needed to know which state
    % is meant with i_X.
    [MinColl,~,ind_traits] = calc_effect_window(par_plot,[filepath,filename{i}],Twin,opt_ecx,opt_conf);
    
    for j = 1:length(MinColl) % run through MFs
        res_tmp = MinColl{j};
        [~,ind_min] = min(res_tmp(:,2)); % find largest effect for this MF
        COLL = cat(1,COLL,[i res_tmp(ind_min,:) ind_traits(ind_min)]);
    end
    waitbar(i/length(filename),f) % update the waiting bar
end
COLL = sortrows(COLL,[1 2]); % make sure it is sorted on file and MF (MF may be in wrong order)
close(f) % close the waiting bar

%% Plot the results in a single multi-panel plot
% This might not be such a brilliant idea if someone wants to rund hundreds
% of profiles.

n = ceil(sqrt(length(filename)));
[figh,ft] = make_fig(n,n,2); % create a figure window of correct size

for i = 1:length(filename)
    hs = subplot(n,n,i);
    set(hs,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
    if n>1 % only shrink white space when there are more than 1 rows and columns
        p = get(hs,'position');
        p([3 4]) = p([3 4])*1.1; % add 10 percent to width and height
        set(hs, 'position', p);
    end
    
    hold on
    title(filename{i},'interpreter','none',ft.name,ft.title)
    ind_i = COLL(:,1)==i ;
    
    area([MF(1) mf_crit],[1 1],'FaceColor','y','LineStyle','none','FaceAlpha',0.3)
    plot([MF(1) MF(end)],[0.9 0.9],'k--')
    
    plot(COLL(ind_i,2),COLL(ind_i,3),'ko:','MarkerFaceColor','r')
    set(hs, 'XScale', 'log')
    if i > length(filename) - n % only x-labels in last row
        xlabel('MF',ft.name,ft.label)
    end
    if (i-1)/n == floor((i-1)/n) % only y-label in first column
        ylabel('smallest respons',ft.name,ft.label)
    end
    ylim([0 1])
    xlim([MF(1) MF(end)])
end

%% Display summary results on screen

% see which states are there; initialise with silly high values so they
% won't be set unless they exist (and switch-case runs well)
locS  = 100;
locL  = 101;
locR  = 102;
loc_h = 103;
if isfield(glo,'locS') % then we have a state of survival
    locS = glo.locS; % collect the location
end
if isfield(glo,'locL') % then we have a state of body length
    locL = glo.locL; % collect the location
end
if isfield(glo,'locR') % then we have a state of reproduction
    locR = glo.locR; % collect the location
end
if isfield(glo,'loc_h') % then we have a state of healthy
    loc_h = glo.loc_h; % collect the locations
end

diary(glo.diary) % collect output in the diary "results.out"
disp(' ')
disp('=================================================================================')
disp('Filename of profile          range for minimum EP10    based on endpoint/window')
disp('=================================================================================')

prof_interest = [];
twin = [];
for i = 1:length(filename)
    fprintf('%-40s ',filename{i})
    COLL_tmp = COLL(COLL(:,1)==i,:);
    ind1 = find(COLL_tmp(:,3)>0.9,1,'last');
    ind2 = find(COLL_tmp(:,3)<=0.9,1,'first');
    tst = 0;
    if ind1 == size(COLL_tmp,1)
        fprintf('  >%4.0f   ',COLL_tmp(end,2))
        ind2 = ind1;
        tst  = 1; % flag to warn that there are only small effects
    elseif ind2 == 1
        fprintf('  <%4.0f   ',COLL_tmp(1,2))
        ind1 = ind2;
        tst  = 2; % flag to warn that there are only large effects
    else
        fprintf('%4.0f-%4.0f ',COLL_tmp([ind1 ind2],2))
    end
    if tst ~= 1 % only if there is sufficient effects
        fprintf('        ')
        switch COLL_tmp(ind2,5)
            case -1
                fprintf('intr. rate')
            case locS
                fprintf('survival')
            case locL
                fprintf('length')
            case locR
                fprintf('repro')
            case loc_h
                fprintf('healthy')
        end
        fprintf(' (t=%1.0f)',COLL_tmp(ind2,4))
    else
        tst = 0; % reset it for next profile
    end
    
    if COLL_tmp(ind1,2) < mf_crit
        fprintf(' *')
        prof_interest = cat(1,prof_interest,COLL_tmp(1,1));
        twin          = cat(1,twin,COLL_tmp(ind2,4));
    end
    
    fprintf('\n')
end
disp('=================================================================================')

% If the profile is not located in the working directory, return the partial path
if ~strcmp(filepath(1:end-1),pwd)
    a = strfind(filepath,filesep);
    if length(a) > 2
        disp(['profiles from directory: ',filepath(a(end-2):end)])
    else
        disp(['profiles from directory: ',filepath])
    end
end

file_out = {};
if ~isempty(prof_interest)
    disp(['* marks profiles that may warrant a closer look (EP10<',num2str(mf_crit),')'])
    disp('  note that the reported window is not necessarily the one with the lowest EP10')
    disp(' ')
    disp(['Profiles with EP10<',num2str(mf_crit),':'])
    disp('========================================================')
    for i = 1:length(prof_interest)
        disp(filename{prof_interest(i)})
        if strcmp(filepath(1:end-1),pwd) % collect the filenames of the interesting profiles
            file_out{i} = filename{prof_interest(i)};
        else % if the profile is not located in the working directory, return the entire path
            file_out{i} = [filepath,filename{prof_interest(i)}];
        end
    end
    disp('========================================================')
end
disp(' ')
diary off

