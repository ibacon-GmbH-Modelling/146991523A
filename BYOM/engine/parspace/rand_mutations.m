function coll_tries = rand_mutations(pmat,coll_ok,bnds_tmp,n_tr_i,d_grid_i,f,n_rnd)

% Usage: coll_tries = rand_mutations(pmat,coll_ok,bnds_tmp,n_tr_i,d_grid_i,f,n_rnd)
% 
% This function randomly mutates the parameter sets in <coll_ok>, with the
% settings obtained for this round from <calc_parspace>. For each mutated
% set, it calculates the min-log-likelihood (MLL). Mutated sets with their
% MLL are returned in the matrix <coll_tries>.
% 
% The inputs may be organised in a different form. Quite a few inputs are
% needed, and it is only called from <calc_parspace>. However, this part of
% the code is called at two locations in <calc_parspace>, and is therefore
% best placed in a separate function.
% 
% Inputs
% <pmat>     parameter matrix
% <coll_ok>  parameter sets (and their MLL) that need to be mutated to new tries (matrix)
% <bnds_tmp> boundaries of parameter space (log-transformed were needed) (matrix)
% <n_tr_i>   number of new parameter tries per element of <coll_ok>
% <d_grid_i> initial grid spacing for each parameter times the max jump size in this round (vector)
%            Note: the function call in <calc_parspace> produces this entry
%            as <f_d_i*d_grid>.
% <f> and <n_rnd> are only used to update the progress bar
% 
% Outputs
% <coll_tries>  mutated parameter sets with their MLL (matrix)
% 
% Author     : Tjalling Jager
% Date       : May 2020
% Web support: <http://www.openguts.info> and <http://www.debtox.info/byom.html>

% =========================================================================
% Copyright (c) 2018-2021, Tjalling Jager (tjalling@debtox.nl). This file
% is a slightly modified version of the <rand_mutations> code that is
% distributed as part of the Matlab version of openGUTS (see
% http://www.openguts.info). Therefore, this code is distributed under
% the same license as openGUTS (GPLv3). The modifications are not in the
% algorithm itself but only to ensure that the code operates in the general
% BYOM framework.
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%  
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% =========================================================================

%% BLOCK 1. Initial things.

n_fit      = sum(pmat(:,2));  % number of fitted parameters
n_cont     = size(coll_ok,1); % how many sets to mutate with
coll_tries = 1./zeros(n_cont*n_tr_i,n_fit+1); % initialise matrix (with INFs) to catch all new sets and their MLL

%% BLOCK 2. Fill coll_tries.
% This section loops over all parameter sets in <coll_ok>. The first
% sub-loop creates a series of mutated parameter sets, and the second
% sub-loop runs it through <transfer> to obtain an MLL.

for i_ok = 1:n_cont % run through the ok parameter sets that were entered
    % Make new tries out of this parameter set, with a random mutation for
    % each parameter, the mutation jumps at max a factor <f_d_i> of the
    % initial grid spacing in round 1 (max jump distance for each parameter
    % in this round is already integrated in the vector <d_grid_i>).
    
    waitbar(i_ok/n_cont,f,['Round ',num2str(n_rnd),' of mutations']) % update progress bar
    
    % BLOCK 2.1. Generate a matrix with mutated parameter sets for this
    % entry in <coll_ok>.
    p_try = zeros(n_tr_i,n_fit); % initialise a fresh matrix for the new parameter tries
    for i_p = 1:n_fit % run through all fitted parameters
        % Create a vector with <n_tr_i> new parameter values to try, and place in correct column of <p_try>.
        p_try(:,i_p) = (coll_ok(i_ok,i_p) + (rand(1,n_tr_i) * 2 - 1) * d_grid_i(i_p))'; % add a random number between -1 and 1, multiplied by max jump
        p_try(:,i_p) = max(p_try(:,i_p),bnds_tmp(i_p,1)); % make sure that all new tries are within min bound
        p_try(:,i_p) = min(p_try(:,i_p),bnds_tmp(i_p,2)); % make sure that all new tries are within max bound
    end
    
    % BLOCK 2.2. Calculate a minus-log-likelihood for the mutated values in
    % <p_try>, and add the parameters and MLL to <coll_tries>.
    for i_t = 1:n_tr_i % run through all new random tries for this entry in <coll_ok>
        pfit    = (p_try(i_t,:))'; % put new tries in vector <pfit> for parameters that are to be fitted
        mll_tst = transfer(pfit,pmat); % calculate the min-log-likelihood for this parameter combination
        coll_tries((i_ok-1)*n_tr_i+i_t,:) = [p_try(i_t,:) mll_tst]; % and collect it with the parameter values in <coll_tries>
    end
end

%% BLOCK 3. Prepare <coll_tries> and sort.

coll_tries = coll_tries(~isinf(coll_tries(:,end)),:); % remove the ones that are INF for the MLL
coll_tries = sortrows(coll_tries,n_fit+1); % sort based on MLL (keep rows/sets together)

