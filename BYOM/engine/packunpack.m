function pout = packunpack(sel,par,pmat)

% Usage: pout = packunpack(sel,par,pmat)
%
% A handy function to allow working with parameter structures when fitting
%
% <sel> is the selection of:
%       1) unpack structure <par> to a matrix <pout>
%       2) pack matrix <pmat> into a structure <pout>
%
% The next series of options are for fitting a spline with paramaters that
% specify the location of the 'nodes'. For this purpose, the global
% <glo.spln> is used.
%
%       3) pack matrix <pmat> into existing structure <par> to <pout>
%       4) unpack structure <par> into a matrix <pout>
%
% Author     : Tjalling Jager 
% Date       : January 2020
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 % field names, set in prelim_checks

switch sel
    
    case 1 % from structure to matrix form
        names   = glo2.names;
        nfields = length(names); % number of fields in par

        pout = ones(nfields,1) * [0 0 -inf +inf 1]; % initialise pout as matrix
        for i = 1:nfields % run though all fields
            li = length(eval(['par.',names{i}]));  % length of parameter i
            pout(i,1:li) = eval(['par.',names{i},';']);
        end
        % this way, we can also define a parameter as a 1 vector, and it
        % will get a no-fit sign and a wide range -inf to +inf.
        
    case 2 % from matrix to structure form
        names   = glo2.names;
        nfields = length(names); % number of fields in par

        for i = 1:nfields % run though all fields
            eval(['pout.',names{i},' = pmat(i,:);']); % create structure par
        end
        
    case 3 % from matrix to spline parameters to be fitted
        glo.spln{1} = pmat(:,1); % put the time for the spline in a global
        pout = par; % copy existing structure and add spline nodes as parameters
        for i = 1:size(pmat,1) % run through time points
            eval(['pout.',glo.spln{2},num2str(i),'=[pmat(',num2str(i),',2) 1 ',num2str(glo.spln{3}(1)),' ',num2str(glo.spln{3}(2)),'];']); 
            % create parameters to be fitted
        end
        
    case 4 % extract the spline parameters from the parameter structure       
        pout(:,1) = glo.spln{1}; % and create a matrix for splining in derivatives
        for i = 1:length(pout) % run through time points
            pout(i,2) = eval(['par.',glo.spln{2},num2str(i),'(1)']);
            % second columns of the matrix will be the 'nodes' for the spline
        end

end