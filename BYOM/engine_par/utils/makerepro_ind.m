function [data,w] = makerepro_ind(R,opt)

% Usage: [data,w] = makerepro_ind(R,opt)
%
% A small function to rework a reproduction data set with data for each
% individual separatley into a format to use in BYOM. For each individual,
% the number of neonates observed must be entered in <R> (not cumulated:
% just like observed in the test), with a -1 where the first eggs are
% observed in the brood pouch as well as for 'true zeros' later in life
% (moults without neonates). First column should be time (in days) and the
% first row should be the treatment identifiers. 
% 
% When a mother dies, the subsequent repro data should be given a NaN
% rather than a zero. There is room for discussion about what to do with
% any eggs or neonates that are observed at the same time at which the
% mother is found dead. I propose to count them as is, and use NaN for the
% next time points. 
% 
% Output depends on the settings in <opt>: 
% 0) Check whether we can use a single intermoult period for the entire
%    data set. Screen output will show mean intermoult times and brood 
%    sizes across the replicates, as function of brood number and treatment.
% 1) Cumulate reproduction, but remove the time points with zero
%    reproduction. Good for clutch-wise reproduction.
% 2) Cumulate reproduction, but don't remove zeros. Good for continuous
%    reproduction.
% 3) Shift neonate release back to the previous moult. When this option is
%    used, don't shift the model predictions with <glo.Tbp>: the data now
%    represent egg production rather than neonate release.
%
% Author     : Tjalling Jager 
% Date       : June 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

Rout    = R;              % make the data set equal to R (but will be changed below)
R       = R(2:end,2:end); % only data, without times and treatments
t       = Rout(2:end,1);  % time vector
all_id  = Rout(1,2:end);  % all id's, with replication
c_u     = unique(all_id); % unique treatment identifiers
        
switch opt
    
    case 0 % calculate the intermoult for each replicate
    
        max_broods = 0;     % this will catch the max number of broods
        for i = 1:size(R,2) % run through individuals
            
            Rtmp     = R(:,i);       % extract one individual
            ind_eggs = find(Rtmp == -1); % index for time with first egg
            ind_moults = (Rtmp~=0 & ~isnan(Rtmp)); % find observations where repro is not zero (could be -1)           
            Rtmp(ind_eggs) = 0;  % turn the -1 into a zero for brood size means
            
            interm{i}  = diff(t(ind_moults)); % intermoult periods
            brdsz{i}   = Rtmp(ind_moults); % collect brood size
            brdsz{i}   = brdsz{i}(2:end); % remove first one (which is a zero for the first egg)
            max_broods = max(max_broods,length(interm{i})); % maximum number of broods across entire data set
        end
        
        interm_out = nan(max_broods,size(R,2)); % initialise matrix with all intermoults
        brdsz_out  = nan(max_broods,size(R,2)); % initialise matrix with all intermoults
        for i = 1:size(R,2) % run through individuals
            interm_out(1:length(interm{i}),i) = interm{i};
            brdsz_out(1:length(brdsz{i}),i)   = brdsz{i};
        end
          
        disp(' ')
        disp('Checking the data set regarding the intermoult period')
        disp(['Overall mean intermoult time  : ',num2str(mean(interm_out(:),'omitnan'))])
        disp(['Overall std of intermoult time: ',num2str(std(interm_out(:),'omitnan'))])
        disp(' ')
        
        Rmim           = nan(max_broods+1,length(c_u)+1); % initialise matrix for mean intermoults
        Rmim(1,:)      = [-1 c_u];      % first row is identifiers
        Rmim(2:end,1)  = 1:max_broods; % first column is brood nr
        Rbrdsz         = Rmim; % same matrix for brood size
        for i = 1:length(c_u)               % run through treatments
            ind_c_u = find(all_id==c_u(i)); % find replicates for treatment i
            Rmim(2:end,i+1) = mean(interm_out(:,ind_c_u),2,'omitnan'); % mean over replicates, for each brood
            Rbrdsz(2:end,i+1) = mean(brdsz_out(:,ind_c_u),2,'omitnan'); % mean over replicates, for each brood
        end
        % Rout is now the mean intermoult duration for each treatment, for
        % each brood. With this matrix, you can see if intermoult time
        % depends on brood number, and if it depends on the treatments.
        
        disp('Mean intermoult period across broods and across treatments')
        disp('------------------------------------------------------------')
        disp(Rmim)
                
        disp('Mean brood size across broods and across treatments')
        disp('------------------------------------------------------------')
        disp(Rbrdsz)
        
        Rout = 0; % clear Rout, so it is not accidentally used for DATA
        
    case 1 % remove zeros and cumulate (-1 should be used for the appearance of first eggs)
        
        warning('off','backtrace') % no need to display where the warning is generated
        for i = 1:size(R,2)          % run through individuals
            Rtmp     = R(:,i);       % extract one individual
            
            ind_dth = find(isnan(Rtmp),1,'first'); % see if the individual died in this replicate
            if isempty(ind_dth) % if it does not die ...
                ind_dth = length(t); % take last time point
            end
                
            % We may want to check the intermoult period: a long intermoult
            % (based on repro events) may indicate that we need to add a
            % true zero reproduction.
            ind_moults = (Rtmp~=0  & ~isnan(Rtmp)); % find observations where repro is not zero (could be -1)           
            interm     = diff(t(ind_moults)); % intermoult periods
            
            ind_zero = (Rtmp==0);    % find the zeros in the repro data and remember them
            ind_eggs = find(Rtmp == -1); % index for time with first egg, or moults without neonates
            
            if isempty(ind_eggs)     % if there is no observation on first egg ...
                ind_eggs = 1;        % assume that only the first observation time is a true zero
            else
                Rtmp(ind_eggs) = 0;  % turn the -1 into a zero for cumsum
            end
            
            cRtmp    = cumsum(Rtmp);  % cumulative reproduction
            cRtmp(ind_zero) = NaN;    % replace zero repro observations by NaNs
            cRtmp(1:ind_eggs(1)) = 0; % point where first egg is observed is true zero
            
            R(:,i)   = cRtmp;        % put the cumulated repro back into R in the correct column
            
            ind_last = find(~isnan(cRtmp)==1,1,'last');
            interm_last = t(ind_dth) - t(ind_last); % time between last repro event and end of test
            % We may want to check this when it is much longer than the
            % general intermoult period: may indicate that we need to add a
            % true zero reproduction.
            
            if ~isempty(interm) && (any(interm > 5) || interm_last > 5)
                warning(['Some true zeros may need to be added for individual ',num2str(i), ', which is in treatment ',num2str(all_id(i))])
                % error
            end
            
        end
        warning('on','backtrace')
        
        Rout(2:end,2:end) = R; % use the new R for the output, to be used in the DATA array
        
    case 2 % just cumulate and don't remove any zeros (generally not advised)
        
        for i = 1:size(R,2)          % run through individuals
            Rtmp     = R(:,i);       % extract one individual
            ind_eggs = find(Rtmp == -1); % index for time with first egg, or moults without neonates
            if ~isempty(ind_eggs)     % if there such observations on first egg/moults ...
                Rtmp(ind_eggs) = 0;   % turn the -1 into a zero for cumsum
            end
            cRtmp    = cumsum(Rtmp); % cumulative reproduction
            R(:,i)   = cRtmp;        % put the cumulated repro back into R in the correct column
        end
        
        Rout(2:end,2:end) = R; % use the new R for the output, to be used in the DATA array

    case 3 % shift neonate numbers back to previous moult
        
        for i = 1:size(R,2)          % run through individuals
            Rtmp     = R(:,i);       % extract one individual
            
            ind_moults = find(Rtmp ~= 0  & ~isnan(Rtmp)); % find all NON-zeros in the repro data (these are the approximate moults)
            ind_zero   = find(Rtmp == 0); % find the zeros in the repro data
            ind_eggs   = find(Rtmp == -1); % index for time with first egg, or moults without neonates
            
            if isempty(ind_eggs)     % if there is no observation on first egg ...
                error('Shifting neonate releases back to previous moult requires observation on first egg (and moults without neonates).')
            elseif Rtmp(ind_moults(1)) ~= -1 % also if the first moult is not a -1
                error('Shifting neonate releases back to previous moult requires observation on first egg')
            end
            
            Rtmp(ind_eggs) = 0;      % turn the -1's into zeros for cumsum
            Rtmp(ind_moults(1:end-1)) = Rtmp(ind_moults(2:end)); % shift the data to previous moult
            
            cRtmp    = cumsum(Rtmp); % cumulative reproduction
            ind_zero = unique([ind_zero;ind_moults(end)]); % add the last moult to the zero, so it will be made NaN
            cRtmp(ind_zero) = NaN;   % replace zero observations by NaNs
            
            % Some of the initial observations are true zeros.
            interm = diff(t(ind_moults)); % intermoult periods
            Tegg   = t(ind_eggs(1));      % time point for first egg
            Tzeros = Tegg - interm(1);    % time till which we need true zeros
            cRtmp(t<=Tzeros) = 0;         % point where first egg is observed is true zero
            % this subtracts the first intermoult period that we can know:
            % the time between the appearance of eggs and the first brood
            
            R(:,i)   = cRtmp;        % put the cumulated repro back into R in the correct column
            
        end
        
         Rout(2:end,2:end) = R; % use the new R for the output, to be used in the DATA array
end
  
data = Rout; % data matrix is Rout
w = ones(size(Rout)); % weight matrix of ones (as we have individuals)
w(isnan(data)) = 0; % if the data is NaN, make the weight zero
w = w(2:end,2:end); % remove the column for times and the row with identifiers

% % Next code is for repro observations at points where female was found
% % dead. Assuming that NaN is only entered for the interval AFTER the
% % mother was found dead, we can try to compensate ... this needs thought
% % and checking!
% w(2:end,:) = (w(2:end,:) + w(1:end-1,:))/2;
% Rout(2:end,2:end) = Rout(2:end,2:end) / w;
% data = Rout; % data matrix is Rout
% % This implies that any repro observed in an interval where the mother dies
% % is counted twice (the weight factor is 0.5 there). This might be
% % reasonable, since there may als be mothers that die in an interval in
% % which they could have been able to produce a clutch.

% % In the next part, I subtract 3 days from the time vector for the repro
% % data; reason is that the eggs are produced some 3 days before the
% % offspring are counted. After this shift, the time vector thus represents
% % egg formation, which is the relevant trait for DEB analysis. For now, I
% % think it is better to shift the model output with glo.Tbp instead, or 
% % shift to previous moult (option 3 above).
% DATA{3}(2:end,1) = DATA{3}(2:end,1) - 3; % extract x days from time vector
% DATA{3}([2 3 4],:) = []; % remove first two time points from DATA (t<0)
% W{3}([1 2 3],:)    = []; % remove first two time points from W (t<0)
