% Folders and files to include for System Identification Experiments
%{
--------------------- Description -----------------------------------------
Script that takes a list of folders and files, and makes sure they are
included in the MATLAB path. Run this file before running any System
Identification experiment script so that all relevant files and folders are
included. 

------------------------- Versions ----------------------------------------
%{
v1 : Suraj R Pawar, 6-9-2020
    - Initialize
v2 : Suraj R Pawar, 7-20-2020
    - Added code to set latex as default interpreter
%}
v3 : Suraj R Pawar, 7-21-2020
    - Added Results folder to inclusion list
%}

list = {'Function Files';
        'Measurement Files';
        'Results';
        '../Common Files'};
    
for i = 1 : length(list)
    addpath(genpath(list{i}));
end

% Set latex to default interpreter
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');