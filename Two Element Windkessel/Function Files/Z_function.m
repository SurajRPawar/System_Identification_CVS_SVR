% -------------------------------------------------------------------------
% Function file for calculation of the elements of the sensitivity
% coefficient matrix. This matrix is based on Yao et al. 2003. It is used
% for ranking parameters on their estimability using available measurements
% during each stage of the cardiac cycle. Refer to PhD Research -> SVR
% Phase 2 -> 8
% Pawar, 2-24-2020
%
% Modified for new two element model, 2-1-2021
% -------------------------------------------------------------------------

function Z = Z_function(theta, x, Qvad, stage)

    % Extract parameters 
        A = theta(1);
        B = theta(2);
        E = theta(3);
        Rsvr = theta(4);
        Cao = theta(5);
        Rv = theta(6);
        Pr = theta(7);
        T = theta(8);
        e = theta(9);
        
    % Extract states
        Vbar = x(1);
        Pao = x(2);
        
    % Calculate Z coefficients
        if stage == 1       % Isovolumic
            
            Z = [-(exp(B*(Vbar - Qvad*T)) - 1)*(e - 1), -A*exp(B*(Vbar - Qvad*T))*(Vbar - Qvad*T)*(e - 1), e*(Vbar - Qvad*T), 0, 0, 0, 0;
                 0, 0, 0, (Pao*T)/(Cao*Rsvr^2), -(T*(Qvad - Pao/Rsvr))/Cao^2, 0, 0;
                 0, 0, 0, 0, 0, 0, 0];
 
                
        elseif stage == 2   % Ejection
            
            Z = [ (E*T*e*(exp(B*Vbar) - 1)*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) - (exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv))) - 1)*(e - 1) - (A*B*T*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(exp(B*Vbar) - 1)*(e - 1)^2)/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)), (A*E*T*Vbar*e*exp(B*Vbar)*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) - A*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(e - 1)*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv) + (A*B*T*Vbar*exp(B*Vbar)*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2))), e*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)) - (E*T*Vbar*e^2)/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) + (A*B*T*Vbar*e*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)), 0, 0, (E*T*e*(E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2))/Rv^2 - (A*B*T*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(e - 1)*(E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2))/Rv^2, 0;
                 -(T*(exp(B*Vbar) - 1)*(e - 1))/(2*Cao*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)), -(A*T*Vbar*exp(B*Vbar)*(e - 1))/(2*Cao*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)), (T*Vbar*e)/(2*Cao*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)), (Pao*T)/(Cao*Rsvr^2), -(T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv - Pao/Rsvr))/Cao^2, -(T*(E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2))/(Cao*Rv^2), 0;
                 -((exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv))) - 1)*(e - 1) - (T*(exp(B*Vbar) - 1)*(e - 1))/(2*Cao*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) - (E*T*e*(exp(B*Vbar) - 1)*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) + (A*B*T*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(exp(B*Vbar) - 1)*(e - 1)^2)/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)))/(2*Rv*(- Pao + E*e*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv)) - (T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv - Pao/Rsvr))/Cao - A*(exp(B*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv))) - 1)*(e - 1))^(1/2)), ((A*T*Vbar*exp(B*Vbar)*(e - 1))/(2*Cao*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) - A*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(e - 1)*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv) + (A*B*T*Vbar*exp(B*Vbar)*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2))) + (A*E*T*Vbar*e*exp(B*Vbar)*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)))/(2*Rv*(- Pao + E*e*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv)) - (T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv - Pao/Rsvr))/Cao - A*(exp(B*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv))) - 1)*(e - 1))^(1/2)), (e*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)) - (T*Vbar*e)/(2*Cao*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) - (E*T*Vbar*e^2)/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)) + (A*B*T*Vbar*e*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(e - 1))/(2*Rv*(- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)))/(2*Rv*(- Pao + E*e*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv)) - (T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv - Pao/Rsvr))/Cao - A*(exp(B*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv))) - 1)*(e - 1))^(1/2)), -(Pao*T)/(2*Cao*Rsvr^2*Rv*(- Pao + E*e*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv)) - (T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv - Pao/Rsvr))/Cao - A*(exp(B*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv))) - 1)*(e - 1))^(1/2)), (T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv - Pao/Rsvr))/(2*Cao^2*Rv*(- Pao + E*e*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv)) - (T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv - Pao/Rsvr))/Cao - A*(exp(B*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv))) - 1)*(e - 1))^(1/2)), ((T*(E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2))/(Cao*Rv^2) + (E*T*e*(E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2))/Rv^2 - (A*B*T*exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)))*(e - 1)*(E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2))/Rv^2)/(2*Rv*(- Pao + E*e*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv)) - (T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv - Pao/Rsvr))/Cao - A*(exp(B*(Vbar - T*(Qvad + (- Pao - A*(exp(B*Vbar) - 1)*(e - 1) + E*Vbar*e)^(1/2)/Rv))) - 1)*(e - 1))^(1/2)) - (E*e*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv)) - Pao - (T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv - Pao/Rsvr))/Cao - A*(exp(B*(Vbar - T*(Qvad + (E*Vbar*e - A*(exp(B*Vbar) - 1)*(e - 1) - Pao)^(1/2)/Rv))) - 1)*(e - 1))^(1/2)/Rv^2, 0];
 
        elseif stage == 3   % Filling
            
            Z = [(E*T*e*(exp(B*Vbar) - 1)*(e - 1))/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)) - (exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv))) - 1)*(e - 1) - (A*B*T*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv)))*(exp(B*Vbar) - 1)*(e - 1)^2)/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)), (A*E*T*Vbar*e*exp(B*Vbar)*(e - 1))/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)) - A*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv)))*(e - 1)*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv) + (A*B*T*Vbar*exp(B*Vbar)*(e - 1))/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2))), e*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv)) - (E*T*Vbar*e^2)/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)) + (A*B*T*Vbar*e*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv)))*(e - 1))/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)),                    0,                            0, (A*B*T*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv)))*(e - 1)*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2))/Rv^2 - (E*T*e*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2))/Rv^2, (E*T*e)/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)) - (A*B*T*exp(B*(Vbar - T*(Qvad - (Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2)/Rv)))*(e - 1))/(2*Rv*(Pr + A*(exp(B*Vbar) - 1)*(e - 1) - E*Vbar*e)^(1/2));
                  0, 0, 0, (Pao*T)/(Cao*Rsvr^2), -(T*(Qvad - Pao/Rsvr))/Cao^2, 0, 0;
                  0, 0, 0, 0, 0, 0, 0];
 
        end
 
end