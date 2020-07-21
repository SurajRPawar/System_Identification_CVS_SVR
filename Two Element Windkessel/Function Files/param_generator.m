function parameters = param_generator(upper_limits, lower_limits)
% Sample parameter values from a uniform distribution
%{
---------------------- Description ----------------------------------------
Sample parameter values of the CVS from a uniform distribution. Values will
be used to generate computational model data and further for system
identification tests. 

--------------------------- Inputs ----------------------------------------
upper_limits : Upper limits for parameters
lower_limits : Lower limits for parameters

-------------------------- Outputs ----------------------------------------
parameters     : Parameters used for computational model data generation
                    (1) A       (mmHg)
                    (2) B       (1/mL)
                    (3) Emax    (mmHg/mL)
                    (4) Cs      (mL/mmHg)
                    (5) Rsvr    (mmHg/mL/s)
                    (6) Rv      (mmHg/mL/s)
                    (7) HR      (bpm)

-------------------------- Versions ---------------------------------------
%{
v1 : Suraj R Pawar, 7-16-2020
    - Initialize
%}
%}

uA = upper_limits(1);
lA = lower_limits(1);
A = lA + (uA - lA)*rand(1);

uB = upper_limits(2);
lB = lower_limits(2);
B = lB + (uB - lB)*rand(1);

uE = upper_limits(3);
lE = lower_limits(3);
E = lE + (uE - lE)*rand(1);

uCs = upper_limits(4);
lCs = lower_limits(4);
Cs = lCs + (uCs - lCs)*rand(1);

uRsvr = upper_limits(5);
lRsvr = lower_limits(5);
Rsvr = lRsvr + (uRsvr - lRsvr)*rand(1);

uRv = upper_limits(6);
lRv = lower_limits(6);
Rv = lRv + (uRv - lRv)*rand(1);

uHR = upper_limits(7);
lHR = lower_limits(7);
HR = floor(lHR + (uHR - lHR)*rand);

parameters = [A; B; E; Cs; Rsvr; Rv; HR];
end