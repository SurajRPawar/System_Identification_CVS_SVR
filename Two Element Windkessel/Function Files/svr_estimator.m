function svr = svr_estimator(noisy_signals)
% Estimated value of Systemic Vascular Resistance
%{
--------------------- Description -----------------------------------------
This function file accepts noisy signals, and calculates the mean Systemic
Vascular Resistance (SVR) using the equations : 
SVR = (Mean Arterial Pressure - Central Venous Pressure)/ (Cardiac Output)

------------------------- Inputs ------------------------------------------
noisy_signals : Signals from the computational model data generator, with
                additive Gaussian white noise added. 
                (1) Plv, Left ventricular pressure (mmHg)
                (2) Pao, Aortic pressure (mmHg)
                (3) Pla, Left Atrial pressure (mmHg)
                (4) Pra, Right Atrial pressure (mmHg)
                (5) Vlv, Left ventricle volume (mL)
                (6) Qa, Aortic flow rate (mL/s)
                (7) Qvad, VAD flow rate (mL/s)
                (8) Time vector (s)

--------------------------- Outputs ---------------------------------------
svr : Calculated value of SVR (mmHg/mL/s)

---------------------------- Versions -------------------------------------
%{
v1 : Suraj R Pawar, 7-15-2020
    - Initialize
%}
%}

Map = mean(noisy_signals(:,2));
Cvp = mean(noisy_signals(:,4));
CO = mean(noisy_signals(:,6) + noisy_signals(:,7));
svr = (Map - Cvp)/CO;
end

