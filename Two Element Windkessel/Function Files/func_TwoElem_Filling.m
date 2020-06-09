function [sigmas_x_new, est_meas] = func_TwoElem_Filling(t, dt, sigmas,...
                                              Qvad, parameters, counts, version)
% Model file to propagate sigma points through filling dynamics
%{
------------------------- DESCRIPTION -------------------------------------
Two Element function file for fillingstage. USE ONLY WITH UKF
Takes the current sigma points, inputs (Qvad) and propagates the sigma
points through filling dynamics to output states and measurements.
This function file uses Euler integration.

----------------------------- INPUTS --------------------------------------
t         : Current time (sec)
dt        : Step size for Euler Integration (sec) 
sigmas    : Current sigma points (include augmented states, process 
                                  and measurement noise terms)
Qvad      : Input signal of VAD flow (mL/s)
parametes : Known parameters
              Rsvr (mmHg.s/mL)
              Cs (mL/mmHg)
              Pr (mmHg)
              Ra (mmHg/mL/s)
              HR (beats per min)
              t_vc : Ventricular contraction time (s)
              t_c  : Total cardiac cycle time (s)
counts    : [num_states; num_meas; num_process_noise_terms; num_meas_noise_terms]

----------------------------- OUTPUTS -------------------------------------
sigmas_x_new : New sigma points propagated through ejection dynamics
est_meas     : Estimated measurements

-------------------------- VERSION HISTORY --------------------------------
v1 : Suraj R Pawar, 5-24-2020
     - Initialize
%}

    % Known parameters
        [getparams, paramdots] = func_handle_parameters(parameters, version, sigmas);
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
        % Normalized elastance
        e_n = normalized_elastance(t, t_c, t_vc);
        
        % Compute pressure using formula from Gohean 2013
        Plv = func_Plv(e_n, A, B, Emax, Vbar);        
        
    % Valve flow
        Qa = 0;
        Qm = (1./Rv).*(Pr - Plv);
        % Do not let Qa be negative
        if any(Qm <= 0) Qm(find(Qm <= 0)) = 0; end
        
    % State Equations for Filling
        Vbardot = -Qvad + Qm;
        Psdot = (1./Cs)*(Qvad - (1./Rsvr)*(Ps - Pr) + Qa);
        
        xdot = [Vbardot; Psdot; paramdots];                
        sigmas_x_new = sigmas_x + xdot*dt + sigmas_v;
                        
        checks;
        
    % Form new measurements
        est_meas(1,:) = Plv + sigmas_n(1,:);
        est_meas(2,:) = sigmas_x_new(2,:) + sigmas_n(2,:);
        est_meas(3,:) = Qa + sigmas_n(3,:);
end
