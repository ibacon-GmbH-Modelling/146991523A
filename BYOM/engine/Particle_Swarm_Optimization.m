%%
% INPUT VARIABLES
% Bird_in_swarm=Number of particle=agents=candidate
% Number_of_quality_in_Bird=Number of Variable
%
% MinMaxRange: jx2 matrix; jth row contains minimum and maximum values of the jth variable
% say you have a variable N1 
% which can have maximum value M1 and minimum value m1
% then your matrix will be [m1 M1] 
% for more:
% [m1 M1; m2 M2; mj Mj]
%
% Food_availability=Objective function with one input variable (for more than one variable you may use array)
% example for two variable
% function f = funfunc(array)
%     a=array(1);
%     b=array(2);
%     f = a+b ;
% end
% Food_availability is a string, for above example : 'funfunc'
%
% availability_type is string 'min' or 'max' to check depending upon need to minimize or maximize the Food_availability
% velocity_clamping_factor (normally 2)
% cognitive_constant=c1=individual learning rate (normally 2)
% social_constant=c2=social parameter (normally 2)
% normally C1+C2>=4
%
% Inertia_weight=At the beginning of the search procedure, diversification is heavily weighted, while intensification is heavily weighted at the end of the search procedure.
% Min_Inertia_weight=min of inertia weight (normally 0.4)
% Max_Inertia_weight=max of inertia weight (normally 0.9)
% max_iteration=how many times readjust the position of the flock/swarm of birds its quest for food
% 
%
% OUTPUT VARIABLE
% optimised_parameters : Optimal parameters
%%

% Downloaded from Matlab Central. Tjalling Jager made small changes to allow
% passing on additional arguments to the objective function, and to modify the output to screen.
% 
% Copyright (c) 2013-14, Pramit Biswas
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [ optimised_parameters ] = Particle_Swarm_Optimization (Bird_in_swarm, Number_of_quality_in_Bird, MinMaxRange, Food_availability, availability_type, velocity_clamping_factor, cognitive_constant, social_constant, Min_Inertia_weight, Max_Inertia_weight, max_iteration,varargin)
disp(' ')
%{
    Checking all functions are present
%}
if (exist ('MinMaxCheck.m')==0)
    % clc;
    fprintf ('Please download the following submission from: <a href="http://www.mathworks.com/matlabcentral/fileexchange/43251-bound-values-of-an-array">MATLAB File Exchange (Click here to open link)</a> \ndownload code by clicking "Download Submission" button \nthen extract and put MinMaxCheck.m in current directory and try again\n');
    return;
end
%{
    Checking all parameteres are entered
%}
if nargin < 11
    error('Missing input parameter(s)!')
end

%{
	universalize availability type  
%}
availability_type=lower(availability_type(1:3));

%{
	 Checking for proper boundary Values and entered Matrix
%}
[row,col]=size(MinMaxRange);
if row~=Number_of_quality_in_Bird || col~=2
    error('Not a proper MinMaxRange Matrix')
end
for i=1:Number_of_quality_in_Bird
    if MinMaxRange(i,1)>=MinMaxRange(i,2)
        error('Minimum value greater than Maximum value!!!')
    end
end

%{
	 counter to display % of completion
%}
N=Bird_in_swarm*max_iteration;
q=0;

%{
	 distinguishing min and max range
%}
bird_min_range=MinMaxRange(:,1);
bird_max_range=MinMaxRange(:,2);

%{
	 
%}
format long;
for i=1:Number_of_quality_in_Bird
    bird(:,i)=bird_min_range(i)+(bird_max_range(i)-bird_min_range(i))*rand(Bird_in_swarm,1);
end

Vmax=bird_max_range*velocity_clamping_factor;
Vmin=-Vmax;

for i=1:Number_of_quality_in_Bird
    Velocity(:,i)=Vmin(i)+(Vmax(i)-Vmin(i))*rand(Bird_in_swarm,1);
end

gBest_availability = +inf;

for itr=1:max_iteration
    % fprintf('Completed  %d %% ...', uint8(q*100/N ))
    % fprintf('Completed %-2.2g%% ...', q*100/N)
    
    fprintf('Completed %-2.2g%%, min f(x): %g\n',q*100/N,gBest_availability)
    
    for p=1:Bird_in_swarm
        parameter=bird(p,:,itr);
        availability(p,itr)=feval(Food_availability,parameter,varargin{1});
        
        switch availability_type
            case 'min'
                format long;
                [pBest_availability,index]=min(availability(p,:));
                pBest=bird(p,:,index);
                
                if(p==1 && itr==1)
                    gBest=pBest;
                    gBest_availability=pBest_availability;
                elseif availability(p,itr)<gBest_availability
                    gBest_availability=availability(p,itr);
                    gBest=bird(p,:,itr);
                end
                
            case 'max'
                format long;
                [pBest_availability,index]=max(availability(p,:));
                pBest=bird(p,:,index);
                
                if(p==1 && itr==1)
                    gBest=pBest;
                    gBest_availability=pBest_availability;
                elseif availability(p,itr)>gBest_availability
                    gBest_availability=availability(p,itr);
                    gBest=bird(p,:,itr);
                end
                
            otherwise
                error('availability_type mismatch')
        end
        
        w(itr)=((max_iteration - itr)*(Max_Inertia_weight - Min_Inertia_weight))/(max_iteration-1) + Min_Inertia_weight;
        Velocity(p,:,(itr+1))=w(itr)*Velocity(p,:,itr) + social_constant*rand(1,Number_of_quality_in_Bird).*(gBest-bird(p,:,itr)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest-bird(p,:,itr));
        Velocity(p,:,(itr+1))=MinMaxCheck(Vmin, Vmax, Velocity(p,:,(itr+1)));
        
        bird(p,:,(itr+1))= bird(p,:,itr) + Velocity(p,:,(itr+1));
        bird(p,:,(itr+1))=MinMaxCheck(bird_min_range, bird_max_range, bird(p,:,(itr+1)));
        q=q+1;
        
    end
    
    % clc;
    % fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
end
optimised_parameters=gBest;