function data = noisy_twoelem_data(filename)
% Loads noisy data from a two element model simulation
%{
--------------------------- DESCRIPTION -----------------------------------
If the filename is empty, the default measurement file
"Noisy_TwoElem_Data.mat" is loaded.

----------------------------- I/O -----------------------------------------
Inputs
filename : name of the file to load. (string)

Outputs
data : strcture containing all loaded variables. See 'Parameters loaded'
       section below for description of parameters loaded.
----------------------------- Extra Notes ---------------------------------
1. Plv = (1 - e_n)*A*(exp(B*(Vlv - V0)) - 1) + e_n*E*(Vlv - V0)
   where, 
   A, B, E : parameters
   en      : Normalized elastance
2. Ventricular contraction time (t_vc) is from beginning of R wave, to
   when normalized elastance becomes 0.
3. States are x = [Vlv; Ps]

--------------------------- Parameters loaded -----------------------------
dt                : Euler integration time step used (s)
parameters  
        1. Cs     : Systemic compliance (mL/mmHg)
        2. Rsvr   : Systemic Vascular Resistance (mmHg/mL/s)
        3. Pr     : Constant pulmonary circulation pressure (mmHg)
        4. Ra     : Resistane of aortic valve (= mitral valve) (mmHg/mL/s)  
        5. Rm     : Resistance of mitral valve (= aortic) (mmHg/mL/s)  
        6. A      : Parameter used for LV pressure generation (mmHg)
        7. B      : Parameter used for LV pressure generation (1/mL)
        8. E      : Maximum Elastance of LV (mmHg/mL)
        9. V0     : Unstressed volume in LV (mL)
        10. HR    : Heart rate (bpm)
        11. t_c   : Time for one cardiac cycle = 60/HR (sec)
        12. t_vc  : Ventricular contraction time (s)
x0                : Initial condition on states
num_beats         : Number of beats included in measurement file
process_noise_std : Process noise standard deviation for each state
meas_noise_std    : Measurement noise standard deviation for each
                    measurement
x                 : States of the two element windkessel model
Plv               : Left ventricular pressure (mmHg)
Ps                : Systemic pressure (mmHg)
Qa                : Aortic valve flow (mL/s)
Qvad              : VAD flow (mL/s)
stage             : Stage of the cardiac cycle. 
                    2 = Ejection
                    3 = Filling
                    1 = Isovolumic expansion / contraction
t                 : Time vector is formed manually

---------------------------- Version History ------------------------------
%{
v1 : 5-22-2020, Suraj R Pawar
     Initialize function file
%}
v2 : Suraj R Pawar, 6-23-2020
    - Function also creates time vector that can be used for simulations
    - Performs a check to make sure data signals and time vector are of
    same length
%} 

    arguments
        filename = 'Noisy_TwoElem_Data_Long.mat';
    end
    
    data = load(filename);   
    
    % Form time vector and make data the same length as time vector    
    num_beats = data.num_beats;
    dt = data.dt;
    tc = data.parameters(11);
    tf = tc*num_beats;
    t = [0:dt:tf];
    
    diff_numel = numel(data.Plv) - numel(t);
    
    if diff_numel > 0
        % More data samples than time stamps
        data.x = data.x(:,[1 : end - diff_numel]);
        data.Plv = data.Plv(:,[1 : end - diff_numel]);
        data.Ps = data.Ps(:,[1 : end - diff_numel]);
        data.Qa = data.Qa(:,[1 : end - diff_numel]);
        data.Qvad = data.Qvad(:,[1 : end - diff_numel]);
        data.stage = data.stage(:,[1 : end - diff_numel]);
    elseif diff_numel < 0
        % More time stamps than data samples
        t = t(1 : end - diff_numel);
    end
    
    data.t = t;
end