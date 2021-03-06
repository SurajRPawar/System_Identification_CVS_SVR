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
ABE version : No numerical version needed. The code uses version 1 and
            version 2 interchangably
version 4 : Estimation of Cs, Pr during all stages of CVS
version 5 : Estimation of only Cs during all stages of CVS
version 6 : Estimation of A and Rv during filling stage
version 7 : Estimation of A, Rv and Pr during filling stage
version 8 : Estimation of B, Emax and Rv during isovolumic contraction and
            ejection
version 9 : Estimation of B, Emax and Pr during ejection
version 10 : Estimation of A, B, Emax and Rv during filling
version 11 : Estimation of Cs and Rv
%}