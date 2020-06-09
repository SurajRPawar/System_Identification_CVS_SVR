function stage = func_switch(Plv, Ps, Pr)
% Determine the stage of the cardiac cycle
%{
----------------------- DESCRIPTION ---------------------------------------
Determines the stage of the cardiac cycle - Ejection / Iso / Filling, based
on the pressures and flows. Stage is determined as follows :
Plv > Pr and Plv > Ps : Ejection
Plv < Pr and Plv < Ps : Filling
Plv > Pr and Plv < Ps : Isovolumic

--------------------------- INPUTS ----------------------------------------
Plv : Left ventricular pressure (mmHg)
Ps : Aortic pressure (mmHg)
Pr : Constant pulmonary circulation pressure (mmHg)

---------------------------- OUTPUTS --------------------------------------
stage : 2 = Ejection
        3 = Filling
        1 = Isovolumic expansion / contraction
        -1 = Erorr (no stage criteria satisfied)

------------------------ VERSION HISTORY ----------------------------------
v1 : Suraj R Pawar, 5-24-2020
     - Initialize
%}

    if (Plv > Pr) && (Plv > Ps)        
        stage = 2;  % Ejection                   
    elseif (Plv < Pr) && (Plv < Ps)
        stage = 3;  % Filling
    elseif (Plv >= Pr) && (Plv <= Ps)
        stage = 1;  % Isovolumic expansion / contraction
    else
        stage = -1; % Error. No criteria satisfied
    end    
end

