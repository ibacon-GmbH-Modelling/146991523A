function [figh,ft] = make_fig(m,n,varargin)

% Usage: [figh,ft] = make_fig(m,n,varargin)
%
% A short little function that creates a figure window of a certain default
% size. This only works when the setting is that figures are NOT docked!
% Figures of size 1x1, 1x2, 2x2, and 2x3 will be plotted such that the
% individual plots are (or at least: should be) equal sized. Larger plots
% receive a large standard size that should be okay in most cases.
%
% The figure handle is returned by the function in figh. It also returns
% font settings for plots, to make sure that all plots have a similar look.
%
% The position setting is [left bottom width height]. You might need to
% adapt these settings to your monitor. Left is the number of pixels from
% the left edge of the screen and bottom from the bottom edge.
%
% Author     : Tjalling Jager
% Date       : May 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Define font properties for plots.
% This is done in this way to make sure that all plots have the same look,
% makes it easy to modify fonts, and allows to adapt the font size for
% small (laptop) screens. 

% ftnm = 'Courier';
ftnm = 'Arial';

ft.name   = {'FontSize','FontName','FontWeight'}; % names of the properties we want to set
ft.label  = {12,ftnm,'Normal'}; % font definition for axis labels
ft.ticks  = {10,ftnm,'Normal'}; % font definition for ticks (values at axes)
ft.legend = {10,ftnm,'Normal'}; % font definition for legends
ft.annot  = {10,'Courier','Normal'}; % font definition for annotations (probably not used anymore)
ft.text   = {10,ftnm,'Normal'}; % font definition for text added (as overall title with file name)
ft.title  = {10,ftnm,'Bold'};   % font definition for panel titles

% use as: <ft.name,ft.ticks> in <set> statements
% Examples:
% set(h_txt,ft.name,ft.text); % use standard formatting for this header
% set(h_leg,ft.name,ft.legend); % use standard font formats for this legend

% On small screens, Matlab does not allow certain figure sizes, so figure
% will be made smaller in make_fig. This also needs smaller fontsizes to
% work well on small screens.
scrsz = get(groot,'ScreenSize'); % get screensize
% Below, a modification that helps when making plots in batch mode (using
% Matlab's parallel toolbox).
if usejava('desktop') == 0 % this setting is zero when run in batch mode (nodisplay)
    scrsz = [1           1        1680        1050]; % just a regular large screen
end 

if scrsz(4) < 870 % this is probably the limit ... (same limit used in <make_fig>)
    ft.label  = {10,ftnm,'Normal'};
    ft.ticks  = {8,ftnm,'Normal'};
    ft.legend = {8,ftnm,'Normal'};
    ft.annot  = {7,'Courier','Normal'};
    ft.text   = {8,ftnm,'Normal'};
    ft.title  = {8,ftnm,'Bold'};
end

%% See if we immediately need to return
if strcmp(get(0,'DefaultFigureWindowStyle'),'docked')
    figh = figure; % new figure for each set of data, with sub-plots for the states ...
    return % when figures are docked, no need to resize (won't work anyway)
end

flag_tktd = 0;
if ~isempty(varargin) 
    if varargin{1} == 1
    figh = -1; % we are then called only to return ft only
    return
    elseif varargin{1} == 2
        flag_tktd = 1; % flag that we need a figure for plot_tktd
    end
end

%% Otherwise make a figure window

% figure size for 1x1, 1x2, 2x2, 2x3, larger square, larger rectangular
figsz = [420 330;932 330;935 749;1455 749;1000 770;1500 770];
figsz2 = [1500 600; 1600 500]; % test for 2 rows but more columns than 3
extrsz = [100 200]; % pixels from top-right screen corner, for top-right corner of figure window

% on small screens, Matlab does not allow certain figure sizes, so just
% make them small enough ...
if scrsz(4) < 870 % this is probably the limit ...
    figsz(3:end,:) = floor(figsz(3:end,:)/1.6);
    figsz2 = floor(figsz2/1.6);
    extrsz = [50 100]; % modify position for small screens
end

figh = figure; % new figure for each set of data, with sub-plots for the states ...

if strcmp(get(0,'DefaultFigureWindowStyle'),'docked') % only if it is docked, we can return now
    return
end

if flag_tktd == 1 % then we need a standard TKTD multiplot
    
    figsz = [min(n*300,1500) min(m*300,770)]; % set a default figure size
    
    % On small screens, Matlab does not allow certain figure sizes, so just make them small enough ...
    if scrsz(4) < 870 % this is probably the limit ...
        figsz  = floor(figsz/1.6); % decrease size
        extrsz = [50 100]; % modify position for small screens
    end
    figh.Position = [scrsz(3)-figsz(1)-extrsz(1) scrsz(4)-figsz(2)-extrsz(2) figsz(1) figsz(2)]; % and give it the new size

else    
    if m == 1
        if n == 1
            figh.Position = [scrsz(3)-figsz(1,1)-extrsz(1) scrsz(4)-figsz(1,2)-extrsz(2) figsz(1,1) figsz(1,2)]; % and give them a default size for 1x1
        elseif n == 2
            figh.Position = [scrsz(3)-figsz(2,1)-extrsz(1) scrsz(4)-figsz(2,2)-extrsz(2) figsz(2,1) figsz(2,2)]; % and give them a default size for 1x2
        end
    elseif m == 2
        if n == 2
            % figh.Position = [380 170 850 682]; % and give them a default size for 2x2
            % figh.Position = [380 170 935 749]; % and give them a default size for 2x2
            figh.Position = [scrsz(3)-figsz(3,1)-extrsz(1) scrsz(4)-figsz(3,2)-extrsz(2) figsz(3,1) figsz(3,2)]; % and give them a default size for 2x2
        elseif n ==3
            % figh.Position = [380 170 1325 682]; % and give them a default size for 2x3
            % figh.Position = [380 170 1455 749]; % and give them a default size for 2x3
            figh.Position = [scrsz(3)-figsz(4,1)-extrsz(1) scrsz(4)-figsz(4,2)-extrsz(2) figsz(4,1) figsz(4,2)]; % and give them a default size for 2x3
            % this needs modifcation as it is determined using the
            % sensitivity plots which are modified in size
        elseif n < 6 % bigger than 3 columns but less than 6 ...
            figh.Position = [scrsz(3)-figsz2(1,1)-extrsz(1) scrsz(4)-figsz2(1,2)-extrsz(2) figsz2(1,1) figsz2(1,2)]; % and give them a default size for 3x4 and larger
        else % otherwise ...
            figh.Position = [scrsz(3)-figsz2(2,1)-extrsz(1) scrsz(4)-figsz2(2,2)-extrsz(2) figsz2(2,1) figsz2(2,2)]; % and give them a default size for 3x4 and larger
        end
    else % it is higher than two sub-plots, so pretty big
        if m == n % if it is square
            % figh.Position = [scrsz(3)-1100 scrsz(4)-870 1000 770]; % and give them a default size for 3x3 and larger
            figh.Position = [scrsz(3)-figsz(5,1)-extrsz(1) scrsz(4)-figsz(5,2)-extrsz(2) figsz(5,1) figsz(5,2)]; % and give them a default size for 3x3 and larger
        else % if it is rectangular
            % figh.Position = [380 170 1500 770]; % and give them a default size for 3x4 and larger
            figh.Position = [scrsz(3)-figsz(6,1)-extrsz(1) scrsz(4)-figsz(6,2)-extrsz(2) figsz(6,1) figsz(6,2)]; % and give them a default size for 3x4 and larger
        end
    end
    
end
