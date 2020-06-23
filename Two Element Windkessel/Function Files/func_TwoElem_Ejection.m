function [sigmas_x_new, est_meas] = func_TwoElem_Ejection(t, dt, sigmas,...
                                              Qvad, parameters, counts, version)
% Model file to propagate sigma points through ejection dynamics
%{
------------------------- DESCRIPTION -------------------------------------
Two Element function file for ejection stage. USE ONLY WITH UKF
Takes the current sigma points, inputs (Qvad) and propagates the sigma
points through Ejection dynamics to output states and measurements.
This function file uses Euler integration.

----------------------------- INPUTS --------------------------------------
t         : Current time (sec)
dt        : Step size for Euler Integration (sec) 
sigmas    : Current sigma points (include augmented states, process 
                                  and measurement noise terms)
Qvad      : Input signal of VAD flow (mL/s)
parametes : All parameter values
             1. Rsvr : Systemic Vascular Resistance (mmHg/mL/s)
             2. Cs   : Systemic compliance (mL/mmHg)
             3. Pr   : Constant pulmonary circulation pressure (mmHg)
             4. Ra   : Aortic valve resistance (mmHg/mL/s)
             5. Rm   : Mitral valve resistance (mmHg/mL/s)
             6. HR   : Heart rate (bpm)
             7. t_vc : Ventricular contraction time (s)
             8. t_c  : Cardiac cycle time (s)
             9. A    : Parameter for PLV generation (mmHg)
             10. B   : Parameter for PLV generation (1/mL)
             11. E   : Maximum LV elastance (mmHg/mL)
             12. V0  : Unstressed blood volume (mL)
counts    : [num_states; num_meas; num_process_noise_terms; num_meas_noise_terms]
version    : The version of experiment to run. This verison will change the
             way parameters and their dynamic equations are set. Refer to
             the MATLAB file 'experiment_versions.m'

----------------------------- OUTPUTS -------------------------------------
sigmas_x_new : New sigma points propagated through ejection dynamics
est_meas     : Estimated measurements

-------------------------- VERSION HISTORY --------------------------------
%{
v1 : Suraj R Pawar, 5-5-2020
     - Initialize
v2 : Suraj R Pawar, 5-24-2020
     - Revise the noise model. Sigma_v is now generated using Brownian motion
       based on Ito's formulation   
v3 : Suraj R Pawar, 5-28-2020
    - Add version mechanism to extract known parameters and set dynamic
    equations for estimated parameters based on the experiment being run
v4 : Suraj R Pawar, 6-11-2020
    - Passing variable 'counts' to func_handle_parameters
%}
v5 : Suraj R Pawar, 6-22-2020
    - Switched to square root law for valve resistance
%}

    % Get known / estimated parameters depending on version
        [getparams, paramdots] = func_handle_parameters(parameters, version, counts, sigmas);
        Rsvr = getparams.Rsvr;
        Cs = getparams.Cs;
        Pr = getparams.Pr;       
        Rv = getparams.Rv;
        HR = getparams.HR;
        t_vc = getparams.tvc;
        t_c = getparams.tc;
        A = getparams.A;
        B = getparams.B;
        Emax = getparams.E;
        
    % Count variables
        num_states = counts(1);                 % States of the model + Parameters to be estimated
        num_meas = counts(2);
        num_process_noise_terms = counts(3);    % Terms for states and parameters to estimate
        
    % Extract Current States
        sigmas_x = sigmas([1:num_states],:);    
        sigmas_v = sigmas([num_states+1:num_states + num_process_noise_terms],:);        
        sigmas_n = sigmas([num_states + num_process_noise_terms+1 :end],:);        
        Vbar = sigmas_x(1,:);                           
        Ps = sigmas_x(2,:);        
    
    % Left Ventricle pressure
        % Normalized Elastance        
        e_n = normalized_elastance(t, t_c, t_vc);
                
        % Compute pressure using formula from Gohean 2013
        Plv = func_Plv(e_n, A, B, Emax, Vbar);        
        
    % Valve flow
        
        dPsys = Plv - Ps;       
        if any(dPsys < 0) dPsys(find(dPsys < 0)) = 0; end        
        Qa = (1./Rv).*sqrt(dPsys);                
        
    % State Equations for Ejection
        Vbardot = -Qvad - Qa;
        Psdot = (1./Cs)*(Qvad - (1./Rsvr)*(Ps - Pr) + Qa);
        
        xdot = [Vbardot; Psdot; paramdots];                
        sigmas_x_new = sigmas_x + xdot*dt + sigmas_v;
        
        checks;
        
    % Form new measurements
        est_meas(1,:) = Plv + sigmas_n(1,:);
        est_meas(2,:) = sigmas_x_new(2,:) + sigmas_n(2,:);
        est_meas(3,:) = Qa + sigmas_n(3,:);
end
