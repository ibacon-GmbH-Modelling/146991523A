function save_plot(fig,savenm,varargin)

% Usage: save_plot(fig,savenm,varargin)
%
% For saving figure windows in various styles. As inputs the figure handle
% <fig>, the name to save the plot with in <savenm>. The options for all
% plots can be set with a series of global setting. Temporary settings can
% be used with varargin.
%
% varargin: 1) handle for the text that is used as title in multiplots
%           2) temporary override of glo.saveplt (to make ad hoc plots in script)
%           3) temporary override of glo.saveplt_str (to make ad hoc plots in script)
%           4) temporary override of glo.saveplt_tight (to make ad hoc plots in script)
%
% Global options (made global as they are used in many functions):
%
% glo.saveplt         : 0) do not save plots 1) save as fig file 2) save as jpeg 3) save as pdf
% glo.saveplt_str     : string (LaTeX) to put title above all plots (works best when titles are suppressed)
% glo.saveplt_ftsiz   : font size for the title in str
% glo.saveplt_tight   : make the plot as tight as possible (especially for multiplots)
% glo.saveplt_notitle : set to 1 to suppress automatic titles above graphs (so you can use str)
%
% Author     : Tjalling Jager
% Date       : April 2019
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo

saveplt = glo.saveplt;
% fill fields that have not been defined yet
if isfield(glo,'saveplt_str')
    str     = glo.saveplt_str; % title above the plot (make sure to set opt_prof.notitle  = 1)
else
    str = []; % title above the plot (make sure to set notitle  = 1)
end
if isfield(glo,'saveplt_ftsiz')
    ftsiz   = glo.saveplt_ftsiz; % font size for the title
else
    ftsiz = 18; % font size for the title
end
if isfield(glo,'saveplt_tight')
    tight   = glo.saveplt_tight; % make the plot as tight as possible (especially for multiplots)
else
    tight = 0; % make the plot as tight as possible (especially for multiplots)
end
if isfield(glo,'saveplt_notitle')
    notitle = glo.saveplt_notitle; % suppress titles
else
    notitle = 0; % set to 1 to suppress title above the plot
end

% use extra optional input argument to override the general settings
if length(varargin)>1
    saveplt = varargin{2};
end
if length(varargin)>2 % if another argument is given, it is a temporary string
    str = varargin{3};
end
if length(varargin)>3 % if another argument is given, it is a temporary tight setting
    tight = varargin{4};
end

if notitle == 1
    figure(gcf)
    title([]) % remove title
    if ~isempty(varargin)
        h_txt = varargin{1}; % handle to text in the title of the figure window
        if iscell(h_txt)
            delete(h_txt{1});
        else
            delete(h_txt);
        end
    end
end 

if tight == 1
    tightfig(fig); % makes the multiplot tighter (removes white around)
end
if ~isempty(str)
    % make a title above the plot (LaTeX format)
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
    text(gca, 0.5, 1,str,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',ftsiz)
end

switch saveplt
    case 1 % save it as a fig file
        saveas(fig,savenm)
    case 2 % save it as a jpeg file
        print(fig,'-djpeg',savenm)
    case 3 % save as PDF
        warning off
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print(fig,savenm,'-dpdf')
        warning on
end
