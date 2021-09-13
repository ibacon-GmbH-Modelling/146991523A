function read_datafiles(fnames,TR,opt)

% A small helper function for DEBtox analyses with the special DEBtox2019
% package. Given a cell array of filenames (each defines a data set), this
% helper function combines them into a single DATA array, W array, and
% defines the table with labels. This makes the code in the main script
% more flexible and more readable. This code requires a certain set up of
% the data files! They need to be functions with the following structure:
% 
% [data,w,LabelTable]= FILENAME(TR,opt,studynr)
% 
% Author     : Tjalling Jager 
% Date       : July 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo DATA W

% start with empty output cell arrays
data_all       = {};
w_all          = {};
LabelTable_all = {};

for i = 1:length(fnames) % run through data files
    eval(['[data,w,LabelTable]=',fnames{i},'(TR,opt,',num2str(i),');']) % read data from file i
    data_all       = cat(1,data_all,data); % add data set to the overal data array
    w_all          = cat(1,w_all,w);       % add weights set to the overal weights array
    LabelTable_all = cat(1,LabelTable_all,LabelTable); % add lables to the overal labels array
end
    
% make data and weights global DATA and W
DATA = data_all; 
W    = w_all;

% create a table with nice custom labels for the legends
glo.LabelTable = LabelTable_all; 

