% Symbolic analysis of two element CVS for sensitivity matrix
%{
----------------------- Description ---------------------------------------
Yao et al. (2003) presents a way to perform parameter identifiability by
obtaining a sensitivity matrix 'Z'. 
Z = partial meas / partial theta. 

We also need to scale the sensitivity matrix according to - 
Z_scaled = (theta est / meas)*(partial meas / partial theta)

Then, the steps for ranking the parameters is as follows - 
1. Calculate the magnitude of each column of Z
2. Select the parameter whose column in Z has the largest magnitude (the
sum of squares of the elements) as the first estimable parameter. 
3. Mark the corresponding column as Xl (l = 1 for the first iteration).
4. Calculate Zlhat, the prediction of the full sensitivity matrix Z, using
the subset of columns Xl : Zlhat = pseudo inverse of Xl * Z
5. Calculate the residual matrix Rl = Z - Zlhat.
6. Calculate the sum of squares of the residuals in each column of RL. 
   The column with the largest magnitude corresponds to the next estimable 
   parameter.
7. Select the corresponding column in Z, and augment the matrix XL by 
   including the new column. Denote the augmented matrix as X_{L +1}.
8. Advance the iteration counter by one and repeat steps 4 to 7 until the 
   column of largest magnitude in the residual matrix is smaller than a 
   prescribed cut-off value.

This script calculates symbolically the expression for the sensitivity
matrix.

--------------------------- Versions --------------------------------------
%{
v1 : Suraj R Pawar, 7-28-2020
%}
v2 : Suraj R Pawar, 2-1-2021
    - Modified for the new (slightly varied) two element windkessel model
%}

clear all; close all; clc;
syms e Qvad Rv Pr Vbar Pao Rsvr Cao A B E T
syms Vbar_kk Pao_kk

% Meaning of each symbol 
%{
e       = Normalized elastance 
Qvad    = VAD flow rate at the kth iteration (mL/s)
Rv      = Valve resistance, taken same for aortic and mitral valves (mmHg.s/mL)
Pr      = Characteristic pulmonary circulation pressure (mmHg)
Vbar    = Vlv - V0 at kth iteration. Vlv is left ventricle volume (state),
          and V0 is the unstressed blood volume (unknown constant) (mL)
Vbar_kk = .... at the k+1 th iteration
Pao     = Aortic pressure at the kth iteration (mmHg)
Pao_kk  = .... at the k+1 th iteration
Rsvr    = Systemic Vascular Resistance (mmHg.s/mL)
Cao     = Aortic compliance (mL/mmHg)
A, B    = Constants for the passive left ventricle elastance curve
E       = Maximum left ventricle elastance (mmHg/mL)
T       = Euler integration time step
%}

Plv  = (1 - e)*A*(exp(B*(Vbar)) - 1) + e*E*(Vbar);
x_kk = [Vbar_kk; Pao_kk];

theta = [A, B, E, Rsvr, Cao, Rv, Pr].';

%% Linearized system matrices

% Ejection
    Vbar_kk = Vbar + (-Qvad - (1/Rv)*sqrt(Plv - Pao))*T;
    Pao_kk = Pao + ((1/Cao)*(Qvad - (1/Rsvr)*(Pao) + (1/Rv)*sqrt(Plv - Pao)))*T;
    
    Plv_kk = (1 - e)*A*(exp(B*(Vbar_kk)) - 1) + e*E*(Vbar_kk);
    Qao_kk = (1/Rv)*sqrt(Plv_kk - Pao_kk);
    
    y = [Plv_kk; Pao_kk; Qao_kk];
    z_ejection = jacobian(y, theta)
    
% Filling
    Vbardot = -Qvad + (1/Rv)*sqrt(Pr - Plv);
    Paodot = (1/Cao)*(Qvad - (1/Rsvr)*(Pao));
    Vbar_kk  = Vbar + (Vbardot)*T;
    Pao_kk   = Pao + (Paodot)*T;

    Plv_kk = (1 - e)*A*(exp(B*(Vbar_kk)) - 1) + e*E*(Vbar_kk);
    Qao_kk = 0;
    
    y = [Plv_kk; Pao_kk; Qao_kk];
    z_filling = jacobian(y, theta)
    
    %{
    f = [Vbardot; Paodot; zeros(6,1)];
    F = jacobian(f, [Vbar; Pao; A; B; E; Cao; Pr; Rv])
    H = jacobian([Plv; Pao], [Vbar; Pao; A; B; E; Cao; Pr; Rv])
    %}
    
% Isovolumic
    Vbar_kk  = Vbar + (-Qvad)*T;
    Pao_kk   = Pao + ((1/Cao)*(Qvad  - (1/Rsvr)*(Pao)))*T;

    Plv_kk = (1 - e)*A*(exp(B*(Vbar_kk)) - 1) + e*E*(Vbar_kk);
    Qao_kk = 0;
    
    y = [Plv_kk; Pao_kk; Qao_kk];
    z_isovolumic = jacobian(y, theta)
    