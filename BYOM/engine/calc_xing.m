function [Xing,Xdir] = calc_xing(a,chicrit)

% A small function to calculate the crossings for a profile likelihood in
% <a> with a criterion <chicrit>. This is placed in a function as it is
% used and in <calc_proflik> and in <calc_likregion> (several times in that
% one).
%
% First column in <Xing> is the x-value for the crossing, and the second a
% flag for when the CI is open on the edge (1 instead of a zero in the
% first or last row in <Xing>). Note that a middle column is generated
% below (direction: -1 for a lower bound of the CI, +1 for an upper) but
% this is removed before returning the output.
%
% Author     : Tjalling Jager
% Date       : June 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

a(:,2) = a(:,2) - chicrit; % rephrase the problem to zero finding

a(a(:,2)==0,2) = -1e-100; % avoid the (theoretical) possibility that the
% sign change is *exactly* caught by the profile

tst    = [0;diff(sign(a(:,2)))]; % where does a change sign?
ind_x  = find(tst~=0); % this is were the zero crossings are

if isempty(ind_x) % if there is no zero crossing at all detected ...
    Xing = [a(1,1) -1 1;a(end,1) 1 1]; % use the outer edges of a
    % and add ones in the last column to signal they are open intervals
else
    Xing   = zeros(length(ind_x),3); % initialise matrix to catch the crossings
    for ix = 1:length(ind_x) % run through all crossings
        Xing(ix,1) = interp1([a(ind_x(ix)-1,2) a(ind_x(ix),2)],[a(ind_x(ix)-1,1) a(ind_x(ix),1)],0,'linear');
        % interpolate between the point just before and just after the sign change
        Xing(ix,2) = sign(tst(ind_x(ix))); % -1 for a lower bound of the CI, +1 for an upper
    end
    if Xing(1,2) == 1 % lowest crossing is an upper boundary ...
        Xing = [[a(1,1) -1 1];Xing]; % add the lowest tested value as lower bound
    end
    if Xing(end,2) == -1 % lowest crossing is an upper boundary ...
        Xing = [Xing;[a(end,1) 1 1]]; % add the highest tested value as upper bound
    end
end

Xdir      = Xing(:,2);
Xing(:,2) = []; % remove the middle column as it is not used anymore