%% Track the version of the experiments
%{
----------------------------- Description ---------------------------------
!! THIS FILE DOES NOT CONTAIN CODE, ONLY CONTAINS EXPERIMENT VERSIONS !!
How to use this file : For every experiment, a main file will be written.
This file will then call the SysID UKF file that will run the loop that
performs UKF. Inside this loop, there will be calls to the function files
for the models of the CVS stages. Inside these stages, parameters need to
be handled. The version number of the experiments will determine how these
files are called. 

This versioning scheme helps handle parameters better. 

------------------------------ Experiment versions ------------------------
version 1 : Estimation of B, Emax during isovolumic contraction and
            ejection
version 2 : Estimation of A during isovolumic expansion + filling stage
version 3 : Estimation of A, B, Emax during iso contraction + ejection
%}