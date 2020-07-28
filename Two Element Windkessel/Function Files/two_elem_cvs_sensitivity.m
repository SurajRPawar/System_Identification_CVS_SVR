% -------------------------------------------------------------------------
% Function file for two element windkessel model with euler integration
% Used for calculation of the sensitivity coefficient matrix. Refer to PhD
% Research -> SVR Phase 2 -> 8
% Pawar, 2-24-2020
% -------------------------------------------------------------------------

function [xdot, y, Z] = two_elem_cvs(t, x, parameters, Qvad)
% Extract Parameters ======================================================

    Cao = parameters(1);
    Rsvr = parameters(2);
    Pr = parameters(3);
    Rv = parameters(4);
    A = parameters(5);
    B = parameters(6);
    E = parameters(7);
    t_c = parameters(8);
    t_vc = parameters(9);
    T = parameters(10);
        
% Extract Current States ================================================== 

    Vbar = x(1);
    Pao = x(2);
    
% Compute left ventricle elastance ========================================

    t_n = mod(t,t_c);   % Normalized time in the cardiac cycle
    if t_n >= 0 && t_n < 2*t_vc/3
        e = 0.5*(1 - cos(3*pi*t_n / (2*t_vc)) );
    elseif t_n >= 2*t_vc/3 && t_n < t_vc
        e = 0.5*( 1 + cos(3*pi*t_n/t_vc - 2*pi) );
    else
        e = 0;
    end
    
% Compute left ventricle pressure =========================================
    Plv = (1 - e)*A*(exp(B*(Vbar)) - 1) + e*E*(Vbar);
    
% Determine stage and set valve flows =====================================

    % Always start with isovolumic contraction
        if t_n == 0
            ejection = 0;
            filling = 0;
            Qa = 0;
            Qm = 0;
        else
            % Aortic valve and ejection
                if Plv > Pao;
                    Qa = (1/Rv)*sqrt(Plv - Pao);
                    % Do not let Qa be negative
                        if Qa <= 0
                            Qa = 0;
                        end
                    ejection = 1;                    % Ejection
                else
                    Qa = 0;    
                    ejection = 0;                    
                end
            % Mitral valve and filling
                if Pr > Plv && ejection ~= 1
                    Qm = (1/Rv)*sqrt(Pr - Plv);
                    filling = 1;                    % Filling
                    % Don't let Qm be negative
                        if Qm <= 0
                            Qm = 0;
                        end
                else
                    Qm  = 0;
                    filling = 0;                    
                end
        end
        
% State Equations =========================================================

    theta = [A; B; E; Rsvr; Cao; Rv; Pr; T; e];
    
    if ejection == 1
        % Ejection
        Vbardot  = -Qvad - Qa;
        Paodot   = (1/Cao)*(Qvad - (1/Rsvr)*(Pao - Pr) + Qa);
        stage = 2;        
        Z = Z_function(theta, x, Qvad, stage);
        %{
        Z = [(E*T*e*(exp(B*Vbar) - 1)*(e - 1))/Rv - (exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv))) - 1)*(e - 1) - (A*B*T*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(exp(B*Vbar) - 1)*(e - 1)^2)/Rv, (A*E*T*Vbar*e*exp(B*Vbar)*(e - 1))/Rv - A*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1)*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv) + (A*B*T*Vbar*exp(B*Vbar)*(e - 1))/Rv),                            e*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)) - (E*T*Vbar*e^2)/Rv + (A*B*T*Vbar*e*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1))/Rv,                               0,                                                                                            0,                                                                                                                                                                                                                                                                                                                                                     (A*B*T*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1)*(Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/Rv^2 - (E*T*e*(Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/Rv^2,                0;
            -(T*(exp(B*Vbar) - 1)*(e - 1))/(Cao*Rv), -(A*T*Vbar*exp(B*Vbar)*(e - 1))/(Cao*Rv), (T*Vbar*e)/(Cao*Rv), (T*(Pao - Pr))/(Cao*Rsvr^2), (T*((Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv - Qvad + (Pao - Pr)/Rsvr))/Cao^2, (T*(Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/(Cao*Rv^2), T/(Cao*Rsvr);
            -((exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv))) - 1)*(e - 1) - (T*(exp(B*Vbar) - 1)*(e - 1))/(Cao*Rv) - (E*T*e*(exp(B*Vbar) - 1)*(e - 1))/Rv + (A*B*T*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(exp(B*Vbar) - 1)*(e - 1)^2)/Rv)/Rv, ((A*T*Vbar*exp(B*Vbar)*(e - 1))/(Cao*Rv) - A*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1)*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv) + (A*B*T*Vbar*exp(B*Vbar)*(e - 1))/Rv) + (A*E*T*Vbar*e*exp(B*Vbar)*(e - 1))/Rv)/Rv, (e*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)) - (T*Vbar*e)/(Cao*Rv) - (E*T*Vbar*e^2)/Rv + (A*B*T*Vbar*e*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1))/Rv)/Rv, -(T*(Pao - Pr))/(Cao*Rsvr^2*Rv), -(T*((Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv - Qvad + (Pao - Pr)/Rsvr))/(Cao^2*Rv), (Pao - (T*((Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv - Qvad + (Pao - Pr)/Rsvr))/Cao + A*(exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv))) - 1)*(e - 1) - E*e*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))/Rv^2 - ((T*(Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/(Cao*Rv^2) + (E*T*e*(Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/Rv^2 - (A*B*T*exp(B*(Vbar - T*(Qvad - (Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1)*(Pao + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/Rv^2)/Rv, -T/(Cao*Rsvr*Rv)];
        %}
        
    elseif filling == 1
        % Filling
        Vbardot  = -Qvad + Qm;
        Paodot   = (1/Cao)*(Qvad - (1/Rsvr)*(Pao - Pr));
        stage = 3;
        Z = Z_function(theta, x, Qvad, stage);
        %{
        Z = [ (E*T*e*(exp(B*Vbar) - 1)*(e - 1))/Rv - (exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv))) - 1)*(e - 1) - (A*B*T*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(exp(B*Vbar) - 1)*(e - 1)^2)/Rv, (A*E*T*Vbar*e*exp(B*Vbar)*(e - 1))/Rv - A*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1)*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv) + (A*B*T*Vbar*exp(B*Vbar)*(e - 1))/Rv), e*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)) - (E*T*Vbar*e^2)/Rv + (A*B*T*Vbar*e*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1))/Rv,                             0,                                       0, (A*B*T*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1)*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/Rv^2 - (E*T*e*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e))/Rv^2, (E*T*e)/Rv - (A*B*T*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)/Rv)))*(e - 1))/Rv;
                    0, 0, 0, (T*(Pao - Pr))/(Cao*Rsvr^2), -(T*(Qvad - (Pao - Pr)/Rsvr))/Cao^2,                                                                                                                                                                                                                           0,                                                                                                        T/(Cao*Rsvr);
                    0, 0, 0, 0, 0, 0, 0];
        %}        
        
    elseif ejection == 0 && filling == 0
        % Isovolumic
        Vbardot  = -Qvad;
        Paodot   = (1/Cao)*(Qvad  - (1/Rsvr)*(Pao - Pr));
        stage = 1;
        Z = Z_function(theta, x, Qvad, stage);
        %{
        Z = [ -(exp(B*(Vbar - Qvad*T)) - 1)*(e - 1), -A*exp(B*(Vbar - Qvad*T))*(Vbar - Qvad*T)*(e - 1), e*(Vbar - Qvad*T), 0, 0, 0, 0;
                0, 0, 0, (T*(Pao - Pr))/(Cao*Rsvr^2), -(T*(Qvad - (Pao - Pr)/Rsvr))/Cao^2, 0, T/(Cao*Rsvr);
                0, 0, 0, 0, 0, 0, 0];
        %}
        
    end
    
    xdot    = [Vbardot; Paodot];
    
% Scaling of sensitivity coefficient matrix ===============================
    theta = [A, B, E, Rsvr, Cao, Rv, Pr];
    if stage == 2
        scaling = [theta./Plv; theta./Pao; theta./Qa];
        Z = Z.*scaling;
    else
        scaling = [theta./Plv; theta./Pao];
        Z([1:2],:) = Z([1:2],:).*scaling;
    end
    
% Outputs =================================================================

    y(1) = Qa;
    y(2) = Qm;
    y(3) = Plv;
    y(4) = e;
    y(5) = Plv - Pao;
    y(6) = stage;
end