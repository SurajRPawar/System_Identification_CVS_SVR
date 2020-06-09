% Function file for two element windkessel model
% Refer to PhD file 11-8-19 (1)

function [xdot, y] = two_elem_cvs(t, x, parameters)
% Extract Parameters ======================================================

    Cs      = parameters(1);
    Rsvr    = parameters(2);
    Pr      = parameters(3);
    Ra      = parameters(4);
    Rm      = parameters(5);
    A       = parameters(6);
    B       = parameters(7);
    E       = parameters(8);
    V0      = parameters(9);
    HR      = parameters(10);
    t_c     = parameters(11);
    t_vc    = parameters(12);    
    
% Extract Current States ================================================== 

    Vlv     = x(1);
    Ps      = x(2);    
    
% Compute left ventricle elastance ========================================

%     t_n     = mod(t,t_c);   % Normalized time in the cardiac cycle
%     if t_n >= 0 && t_n < 2*t_vc/3
%         e_n     = 0.5*(1 - cos(3*pi*t_n / (2*t_vc)) );
%     elseif t_n >= 2*t_vc/3 && t_n < t_vc
%         e_n     = 0.5*( 1 + cos(3*pi*t_n/t_vc - 2*pi) );
%     else
%         e_n     = 0;
%     end
    shift = 2/3;
    e_n = normalized_elastance(t, t_c, t_vc, shift);
% Interpolate QVAD signal =================================================
    
    %Qvad = interp1(t_sim, Z_Qvad, t);
    if t_n>= 0.24 && t_n < 0.54;
        Qvad = (2*35/0.3)*(0.5 - 0.5*cos(2*pi*(t_n-0.24)/0.3));
    else 
        Qvad = 0;
    end;
    
% Compute left ventricle pressure =========================================
    Plv     = (1 - e_n)*A*(exp(B*(Vlv - V0)) - 1) + e_n*E*(Vlv - V0);
    
% Determine stage and set valve flows =====================================

    % Always start with isovolumic contraction
        if t_n == 0
            ejection = 0;
            filling = 0;
            Qa = 0;
            Qm = 0;            
        else
            % Aortic valve and ejection
                if Plv > Ps;
                    Qa = (1/Ra)*(Plv - Ps);
                    % Do not let Qa be negative
                        if Qa <= 0
                            Qa = 0;
                        end
                    ejection = 1; % Ejection                    
                else
                    Qa = 0;    
                    ejection = 0;                    
                end
            % Mitral valve and filling
                if Pr > Plv && ejection ~= 1
                    Qm          = (1/Rm)*(Pr - Plv);
                    filling     = 1;                    % Filling
                    % Don't let Qm be negative
                        if Qm <= 0
                            Qm = 0;
                        end
                else
                    Qm  = 0;
                    filling     = 0;                    
                end
        end
        
% State Equations =========================================================

    if ejection == 1
        % Ejection
        Vlvdot  = -Qvad - Qa;
        Psdot   = (1/Cs)*(Qvad - (1/Rsvr)*(Ps - Pr) + Qa);
        stage = 2;
        
    elseif filling == 1
        % Filling
        Vlvdot  = -Qvad + Qm;
        Psdot   = (1/Cs)*(Qvad - (1/Rsvr)*(Ps - Pr));
        stage = 3;
        
    elseif ejection == 0 && filling == 0
        % Isovolumic
        Vlvdot  = -Qvad;
        Psdot   = (1/Cs)*(Qvad  - (1/Rsvr)*(Ps - Pr));
        stage = 1;
    end
    
    xdot    = [Vlvdot; Psdot];
    
% Outputs =================================================================

    y(1) = Qa;
    y(2) = Qm;
    y(3) = Plv;
    y(4) = e_n;
    y(5) = Plv - Ps;
    y(6) = Qvad;
    y(7) = stage;
end