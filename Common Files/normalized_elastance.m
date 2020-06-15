function en = normalized_elastance(t,t_c,t_vc)
% Calculate normalized elastance of the left ventricle
%{
-------------------------- Description ------------------------------------
Equation used from Gohean 2013 paper. 

---------------------------- Inputs ---------------------------------------
t     : Current time (s)
t_c   : Time for one cardiac cycle (s)
t_vc  : Time for ventricular contraction (s)

----------------------------- Outputs -------------------------------------
en : Normalized Elastance

------------------------------ Versions -----------------------------------
v1 : Suraj R Pawar 5-24-2020
    - Initialize
v2 : Suraj R Pawar, 5-28-2020
    - Add description and comments
v3 : Suraj R Pawar, 5-29-2020
    - Shift no longer optional argument, in fact not an argument
    - Fixed en expressions to get correct normalized elastance curve. With
    previous expressions that used 'shift' variable, there was
    discountinuity in the normalized elastance, perhaps due to rounding off
    errors of shift = 2/3
%}    
    % Normalized time
    tn = mod(t,t_c); 
    
    % Normalized Elastace
    if tn >= 0 && tn < 2*t_vc/3
        en = 0.5*(1 - cos(3*pi*tn / (2*t_vc)) );
    elseif tn >= 2*t_vc/3 && tn < t_vc
        en = 0.5*( 1 + cos((3*pi*tn/t_vc) - 2*pi) );
    else
        en = 0;
    end    
end