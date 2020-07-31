function [sigmas_x_new, est_meas] = func_TwoElem_CsRv(t, dt, sigmas,...
                                              u, parameters, counts, version)
% Model file for time propagation of sigma pts for Cs & Rv estimation
%{
------------------------- DESCRIPTION -------------------------------------
!! USE ONLY WITH UKF CODE !!
Dynamics for Pao for the Two Element Windkessel Model. This file is used
for Cs and Rv estimation.

----------------------------- INPUTS --------------------------------------
t         : Current time (sec)
dt        : Step size for Euler Integration (sec) 
sigmas    : Current sigma points (include augmented states, process 
                                  and measurement noise terms)
u         : Input signals
            1. Plv
            2. Qa
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
sigmas_x_new : New sigma points propagated through CVS dynamics
est_meas     : Estimated measurements

-------------------------- VERSION HISTORY --------------------------------
%{
v1 : Suraj R Pawar, 6-22-2020
     - Initialize
%}
v2 : Suraj R Pawar, 7-31-2020
    - Fixed description. 
    - Instead of u in the inputs, the description mentioned Qao. Replaced
    with the right description
%}

    % Get known / estimated parameters depending on version
        [getparams, paramdots] = func_handle_parameters(parameters, version, counts, sigmas);
        Rsvr = getparams.Rsvr;
        Cs = getparams.Cs;
        Pr = getparams.Pr;               
        Rv = getparams.Rv;
        
    % Count variables
        num_states = counts(1);                 % States of the model + Parameters to be estimated
        num_meas = counts(2);
        num_process_noise_terms = counts(3);    % Terms for states and parameters to estimate
        
    % Extract Current States
        sigmas_x = sigmas([1:num_states],:);    
        sigmas_v = sigmas([num_states+1:num_states + num_process_noise_terms],:);        
        sigmas_n = sigmas([num_states + num_process_noise_terms+1 :end],:);                                          
        Ps = sigmas_x(1,:);             
        
    % Flows
        Plv = u(1);
        Qvad = u(2);
        dPsys = Plv - Ps;       
        if any(dPsys < 0) dPsys(find(dPsys < 0)) = 0; end        
        Qa = (1./Rv).*sqrt(dPsys);   
        Qao = Qa + Qvad;
        
    % State Equations for Ejection        
        Psdot = (1./Cs).*(Qao - (1./Rsvr)*(Ps - Pr));
        
        xdot = [Psdot; paramdots];                
        sigmas_x_new = sigmas_x + xdot*dt + sigmas_v;        
        
    % Form new measurements        
        est_meas(1,:) = sigmas_x_new(1,:) + sigmas_n(1,:);
        est_meas(2,:) = Qa + sigmas_n(2,:);
end
