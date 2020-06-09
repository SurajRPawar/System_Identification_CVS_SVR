% Folders and files to include for System Identification Experiments
%{
--------------------- Description -----------------------------------------
Script that takes a list of folders and files, and makes sure they are
included in the MATLAB path. Run this file before running any System
Identification experiment script so that all relevant files and folders are
included. 

------------------------- Versions ----------------------------------------
v1 : Suraj R Pawar, 6-9-2020
    - Initialize
%}

list = {'Function Files';
        'Measurement Files'
        '../Common Files'};
    
for i = 1 : length(list)
    addpath(genpath(list{i}));
end