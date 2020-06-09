function y = func_update_measurements(t, x, est_meas, parameters, stage, version)
% Calculate updated measurements for UKF Sys-ID of two element model
%{
------------------------ DESCRIPTION --------------------------------------
!! Use only for Sys ID with UKF code for Two Element Windkessel model !!
Takes the new mean as calculated through the measurement update step of the
UKF, and updates the measurements of Plv, Ps and Qa

---------------------- INPUTS ---------------------------------------------
t         : Current time (s)
x         : States [Vbar; Ps; A; B; Emax]    
            Vbar : Vlv - V0 (mL)
            Ps   : (mmHg)
            A    : (mmHg)
            B    : (1/mL)
            Emax : (mmHg/mL)
est_meas  : Estimated value of measurements
            Plvhat : (mmHg)
            Pshat : (mmHg)
            Qahat : (mL/s)
parametes : Known parameters
              Rsvr (mmHg.s/mL)
              Cs (mL/mmHg)
              Pr (mmHg)
              Ra (mmHg/mL/s)
              HR (beats per min)
              t_vc : Ventricular contraction time (s)
              t_c  : Total cardiac cycle time (s)
stage     : Stage of cardiac cycle
            2 = Ejection
            3 = Filling
            1 = Isovolumic expansion / contraction
           -1 = Error / Faulty stage
version    : Version number of the experiment

---------------------------- OUTPUT ---------------------------------------
y : Updated measurements [Plv; Ps; Qa]
    Plv : (mmHg)
    Ps  : (mmHg)
    Qa  : (mL/s)

--------------------------- VERSION HISTORY -------------------------------
v1 : Suraj R Pawar, 5-24-2020
     - Initialize
%}

    % Extract states and parameters
    Vbar = x(1);
    Ps = x(2);
    [getparams, ~] = func_handle_parameters(parameters, version, [], x);
    A = getparams.A;
    B = getparams.B;
    Emax = getparams.E;
    
    % Normalized elastance
    t_c = getparams.tc;
    t_vc = getparams.tvc;    
    e_n = normalized_elastance(t, t_c, t_vc);
    
    % Updated measurements
    Plv = func_Plv(e_n, A, B, Emax, Vbar);
    
    switch stage
        case 2
            % Ejection
            Rv = getparams.Rv;
            Qa = (1/Rv)*(Plv - Ps);
            if Qa < 0 Qa == 0; end            
            y = [Plv; Ps; Qa];
        case 1
            % Isovolumic 
            Qa = est_meas(3);
            y = [Plv; Ps; Qa];
        case 3
            % Filling
            Qa = est_meas(3);
            y = [Plv; Ps; Qa];
        otherwise
            Qa = 0;
    end        
end

