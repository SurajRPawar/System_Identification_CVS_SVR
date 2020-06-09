function Plv = func_Plv(e_n, A, B, Emax, Vbar)
% Calculate left ventricle pressure
%{
----------------------------- Description ---------------------------------
Plv equation based on Gohean 2013.
Plv = (1 -en)*A*(exp(B*(Vlv - V0)) -1) + en*Emax*(Vlv-V0)

------------------------------ Inputs -------------------------------------
e_n  : Normalized elastance
A    : Parameter (mmHg)
B    : Parameter (1/mL)
Emax : Maximum left ventricle Elastance (mmHg/mL)
Vbar : Vlv - V0 (mL), where
       Vlv = Left ventricle volume (mL)
       V0 = Unstressed blood volume in LV (mL)

------------------------------- Outputs -----------------------------------
Plv : Left Ventricle Pressure (mmHg)

----------------------------- Versions ------------------------------------
v1 : Suraj R Pawar, 5-24-2020
    - Initialize
v2 : Suraj R Pawar, 5-28-2020
    - Add documentation
%}
    Plv = (1 - e_n)*A.*(exp(B.*Vbar) - 1) + e_n*Emax.*Vbar;    
end

