function [sigmas_x_new, est_meas] = func_TwoElem_Iso(t, dt, sigmas,...
                                              Qvad, parameters, counts)
% Model file to propagate sigma points through isovolumic dynamics
%{
------------------------- DESCRIPTION -------------------------------------
Two Element function file for isovolumic stage. USE ONLY WITH UKF
Takes the current sigma points, inputs (Qvad) and propagates the sigma
points through Isovolumic dynamics to output states and measurements.
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
        Rsvr = parameters(1);
        Cs = parameters(2);
        Pr = parameters(3);       
        Rv = parameters(4);
        HR = parameters(5);
        t_vc = parameters(6);
        t_c = parameters(7);
        A = parameters(8);
        
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
%         A = sigmas_x(3,:);
        B = sigmas_x(3,:);
        Emax = sigmas_x(4,:);
    
    % Left Ventricle pressure
        % Normalized Elastance        
%         t_n = mod(t,t_c);   % Normalized time in the cardiac cycle
%         if t_n >= 0 && t_n < 2*t_vc/3
%             e_n = 0.5*(1 - cos(3*pi*t_n / (2*t_vc)) );
%         elseif t_n >= 2*t_vc/3 && t_n < t_vc
%             e_n = 0.5*( 1 + cos(3*pi*t_n/t_vc - 2*pi) );
%         else
%             e_n = 0;
%         end
        e_n = normalized_elastance(t, t_c, t_vc);
        
        % Compute pressure using formula from Gohean 2013
        Plv = func_Plv(e_n, A, B, Emax, Vbar);
        %Plv = (1 - e_n)*A.*(exp(B.*Vbar) - 1) + e_n*Emax.*Vbar;
        %if Plv <= 0 Plv(find(Plv <= 0)) = 0; end;
        
    % Valve flow
        Qa = zeros(size(Plv));        
        
    % State Equations for Isovolumic stage
        Vbardot = -Qvad - Qa;
        Psdot = (1/Cs)*(Qvad - (1/Rsvr)*(Ps - Pr) + Qa);
%         Adot = zeros(size(Vbardot));
        Bdot = zeros(size(Vbardot));
        Emaxdot = zeros(size(Vbardot));
        xdot = [Vbardot; Psdot; Bdot; Emaxdot];
%         xdot = [Vbardot; Psdot; Adot; Bdot; Emaxdot];
                
        sigmas_x_new = sigmas_x + xdot*dt + sigmas_v;
        
        % A, B, Emax cannot be negative
        checks;
        
    % Form new measurements
        est_meas(1,:) = Plv + sigmas_n(1,:);
        est_meas(2,:) = sigmas_x_new(2,:) + sigmas_n(2,:);
        est_meas(3,:) = Qa + sigmas_n(3,:);
end
