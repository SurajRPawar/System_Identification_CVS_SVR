function [xdot,y] = MCL_CVS_Gohean_ODE(t,x,rho, parameters, vad)
%% CVS ODE function file
% CVS states: x = [VLA, VLV, PSA, QSA, PST, PSV, VRA, VRV, PPA, QPA, PPT, PPV];

VLA  = x(1);    
VLV  = x(2);    
P_sa = x(3);  
Q_sa = x(4);
P_st = x(5); 
P_sv = x(6);    
VRA  = x(7);        
VRV  = x(8);        
P_pa = x(9);  
Q_pa = x(10);
P_pt = x(11); 
P_pv = x(12);

%% Previous Inertial Model Parameters
HR = parameters(7);
Tc = 60/HR; % cardiac cycle period, sec/beat

% Left Ventricle
ELVmax   = parameters(3); % maximum LV elastance, mmHg/mL
TLV   = (550 - 1.75*HR)/1000; % LV contraction time, s
delta = 0.666;    % elastance peak percent shift
ALV   = parameters(1); % passive elastance coefficient, mmHg
BLV   = parameters(2); % passive elastance coefficient, 1/mL

% left atria
ELAmax = 0.25; % mmHg/mL
TLA    = 0.09; % LA contraction time, s
ALA    = 0.45; % passive elastance coefficient, mmHg
BLA    = 0.05; % passive elastance coefficient, 1/mL

% Right Ventricle
ERVmax = 0.75; % maximum RV elastance, mmHg/mL
TRV = TLV;
ARV = 0.04; % ARA in paper, but this may be a typo b/c of ARA in Right Atria
BRV = 0.05; 

% Right Atria
ERAmax = 0.25;
ARA    = 0.45;
BRA    = 0.05;
TRA    = 0.09;


% systemic arterial system
RSA = 0.15;
CSA = parameters(4);
LSA = 0.0022;
RST = 0.8;    % Arterial Tree Resistance
CST = 1.5;      % Arterial Tree Compliance
RSV = parameters(5) - RSA - RST; % Venous Resistance
CSV = 20; 

% CSA = 1.25; % systemic artery compliance, mL/mmHg
% CAT = 1;    % artery tree compliance, mL/mmHg
% RSA = 0.15; % systemic artery resistance, mmHg s/mL
% RAT = 0.8;  % systemic artery tree resistance, mmHg s/mL

% Pulmonary Circulation
RPA = 0.07;
CPA = 7.5;
LPA = 0.0018; % Pulmonary Artery Inertia
RPT = 0.04;
CPT = 0.5; % Labeled as CSP in paper; Arterial Tree Compliance
RPV = 0.003;
CPV = 20;

% Cannulas
Ro = 0.0054;
Lo = 0.0122;

% % inputs
% PRS = 20; % right side pressure, mmHg
% PVS = 10; % venous system pressure, mmHg

% CVS Initial conditions for zero pressure
VLVo = 5;   % zero pressure volume, mL
VLAo = 10;  % zero pressure volume, mL
VRAo = 10;
VRVo = 5;

% Valve Resistance
Rv = parameters(6); 

% % valve flows
% RMI  = 0.0025; % mitral valve resistance, mmHg s/mL
% RAO  = 0.0025; % aortic valve resistance, mmHg s/mL
% RRS  = 0.05; % right side resistance, mmHg s/mL
L_v = 0.0000003;

% %% Pressures
% et = mod(t,Tc); % time within cardiac cycle
% % eLV
% if et<=(TLV*delta)
%     eLV = 0.5*(1-cos(pi*et/(TLV*delta)));
% elseif (et>(TLV*delta))&&(et<=TLV)
%     eLV = 0.5*(1+cos(pi*(et-TLV*delta)/(TLV*(1-delta))));
% else
%     eLV = 0;
% end
% % eLA
% if (et>=(Tc-TLA))&&(et<=Tc)
%     eLA = 0.5*(1-cos(2*pi*(et-Tc+TLA)/TLA));
% else
%     eLA = 0;
% end
% 
% % eRV
% if et<=(TRV*delta)
%     eRV = 0.5*(1-cos(pi*et/(TRV*delta)));
% elseif (et>(TRV*delta))&&(et<=TRV)
%     eRV = 0.5*(1+cos(pi*(et-TRV*delta)/(TRV*(1-delta))));
% else
%     eRV = 0;
% end
% % eRA
% if (et>=(Tc-TRA))&&(et<=Tc)
%     eRA = 0.5*(1-cos(2*pi*(et-Tc+TRA)/TRA));
% else
%     eRA = 0;
% end
% 
% % Pressures
% PLA = (1-eLA)*ALA*exp(BLA*(VLA-VLAo)-1)+ELAmax*eLA*(VLA-VLAo); %mmHg
% PLV = (1-eLV)*ALV*exp(BLV*(VLV-VLVo)-1)+ELVmax*eLV*(VLV-VLVo);
% PRA = (1-eRA)*ARA*exp(BRA*(VRA-VRAo)-1)+ERAmax*eRA*(VRA-VRAo);
% PRV = (1-eRV)*ARV*exp(BRV*(VRV-VRVo)-1)+ERVmax*eRV*(VRV-VRVo);
% 
% %% Flows
% % Mitral Valve Flow
% QMIdot = (PLA - PLV - Rv^2*QMI^2)/L_v;
% if QMI <= 0
%     QMI = 0;
%     if PLA < PLV
%         QMIdot = 0;
%     end
% end
% 
% % Aortic valve flow
% QAOdot = (PLV - PSA - Rv^2*QAO^2)/L_v;
% if QAO < 0
%     QAO = 0;
%     if PLV < PSA
%         QAOdot = 0;
%     end
% end  
% 
% % Tricuspid Valve Flow
% QTValvedot = (PRA - PRV - Rv^2*QTValve^2)/L_v;
% if QTValve < 0
%     QTValve = 0;
%     if PRA < PRV
%         QTValvedot = 0;
%     end
% end
%     
% % Pulmonary Valve Flow
% QPValvedot = (PRV - PPA - Rv^2*QPValve^2)/L_v;
% if QPValve < 0
%     QPValve = 0;
%     if PRV < PPA
%         QPValvedot = 0;
%     end
% end
% 
% % Systemic Flow
% QST = (PST - PSV)/RST;
% QSV = (PSV - PRA)/RSV;
% 
% % Pulmonic Flow
% QPT = (PPT - PPV)/RPT;
% QPV = (PPV - PLA)/RPV;
% 
% %% ODEs
% % Pressures
% PSAdot = (QAO + Qvad - QSA)/CSA;
% PSTdot = (QSA - QST)/CST;
% PSVdot = (QST - QSV)/CSV;
% PPAdot = (QPValve - QPA)/CPA;
% PPTdot = (QPA - QPT)/CPT;
% PPVdot = (QPT - QPV)/CPV;
% 
% % Flows
% QSAdot = (PSA - PST - RSA*QSA)/LSA;
% QPAdot = (PPA - PPT - RPA*QPA)/LPA;
% 
% % Volumes
% VLAdot = QPV-QMI; 
% VLVdot = QMI-QAO-Qvad;
% % VSAdot = QAO-QSA+Qvad;
% % VSTdot = QSA-QST;
% % VSVdot = QST-QSV;
% 
% VRAdot = QSV-QTValve; 
% VRVdot = QTValve-QPValve;
% % VPAdot = QPV-QPA;
% % VPTdot = QPA-QPT;
% % VPVdot = QPT-QPV;

% From Jeff's model:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Systemic Circulation equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RR = Tc;
tn = mod(t,RR);
shift = delta;
R_m = Rv; R_a = Rv; R_t = Rv; R_p = Rv;

% VAD Flow (Counter pulsatile)
if tn>= 0.24 && tn < 0.54;
%     Qvad = 0;
    Qvad = vad*(2*35/0.3)*(0.5 - 0.5*cos(2*pi*(tn-0.24)/0.3));
else 
    Qvad = 0;
end;

% Left atria pressure
    % Normalized elastance (min 0, max 1) 
        if tn >= RR-TLA && tn <= RR
            e_t_la = 1/2 * (1 - cos(2*pi*(tn-(RR-TLA))/TLA));
        else
            e_t_la = 0;
        end
    % Passive elastance pressure
        P_la_passive = (1-e_t_la) * ALA * (exp(BLA*(VLA-VLAo)) - 1);
    % Active elastance pressure
        P_la_active =  + (e_t_la)*ELAmax*(VLA-VLAo);
    % Left atria pressure
        P_la = P_la_passive + P_la_active;
    
% Left ventricular pressure
    % Normalized elastance (min 0, max 1)
        if tn <= TLV*shift
            e_t_lv = 1/2 * (1 - cos(pi*tn/(TLV*shift)));
        elseif tn > TLV*shift && tn <= TLV
            e_t_lv = 1/2 * (1 + cos(pi*(tn-TLV * shift)/(TLV*(1-shift))));
        else
            e_t_lv = 0;
        end
    % Passive elastance pressure
        P_lv_passive = (1-e_t_lv) * ALV * (exp(BLV*(VLV - VLVo)) - 1);
    % Active elastance pressure
        P_lv_active = (e_t_lv)*ELVmax*(VLV-VLVo);
    % Left ventricular pressure
        P_lv = P_lv_active + P_lv_passive;


% % Mitral valve flow   
%     Q_m_dot = (P_la - P_lv - R_m^2*Q_m^2)/L_v;
%     if Q_m <= 0
%         Q_m = 0;
%         if P_la < P_lv
%             Q_m_dot = 0;
%         end
%     end

% Mitral Valve Flow
if P_lv >= P_la
    Q_m = 0;
else
    Q_m = sqrt(P_la - P_lv)/Rv;
end

% % Aortic valve flow
%     Q_a_dot = (P_lv - P_sa - R_a^2*Q_a^2)/L_v;
%     if Q_a < 0
%         Q_a = 0;
%         if P_lv < P_sa
%             Q_a_dot = 0;
%         end
%     end  

% Aortic valve flow
if P_sa >= P_lv
    Q_a = 0;
else
    Q_a = sqrt(P_lv - P_sa)/Rv;
end  
         
    % Pulmonary vein flow    
        Q_pv = (P_pv - P_la) / RPV; %right side return, flow into LA
    % Left atrial volume ODE
        V_la_dot = Q_pv - Q_m;
    % Left ventricular volume ODE
        V_lv_dot = Q_m - Q_a - Qvad; 
        
% Systemic artery pressure and flow ODEs
    P_sa_dot = (Q_a + Qvad - Q_sa)/CSA;
    Q_sa_dot = (P_sa - P_st - Q_sa*RSA)/LSA;
    
    % Systemic artery tree flow and pressure ODE
    Q_st = (P_st - P_sv) / RST;
    P_st_dot = (Q_sa  - Q_st) / CST;   
    
% % Coronary flow
%     if tn > T_v %diastole
%         R_cor = 37;
%     else %systole
%         R_cor = 80;
%     end
%     Q_cor = (P_sa - P_sv) / R_cor;
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pulmonary Circulation equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Right atria pressure
    % Normalized elastance (min 0, max 1) 
        e_t_ra = e_t_la;
    % Passive elastance pressure
        P_ra_passive = (1-e_t_ra) * ARA * (exp(BRA*(VRA-VRAo)) - 1);
    % Active elastance pressure
        P_ra_active =  + (e_t_ra)*ERAmax*(VRA-VRAo);
    % Left atria pressure
        P_ra = P_ra_passive + P_ra_active;

% Right ventricle pressure
    % Normalized elastance (min 0, max 1) 
        e_t_rv = e_t_lv;
    % Passive elastance pressure
        P_rv_passive = (1-e_t_rv) * ARV * (exp(BRV*(VRV-VRVo)) - 1);
    % Active elastance pressure
        P_rv_active =  + (e_t_rv)*ERVmax*(VRV-VRVo);
    % Left atria pressure
        P_rv = P_rv_passive + P_rv_active;
        
% % Tricuspid valve flow
%     Q_t_dot = (P_ra - P_rv - R_t^2*Q_t^2)/L_v;
%     if Q_t <= 0
%         Q_t = 0;
%         if P_ra < P_rv
%             Q_t_dot = 0;
%         end
%     end
    
% Tricuspid Valve Flow
if P_rv >= P_ra
    Q_t = 0;
else
    Q_t = sqrt(P_ra - P_rv)/Rv;
end  


% % Pulmonary valve flow   
%     Q_p_dot = (P_rv - P_pa - R_p^2*Q_p^2)/L_v;
%     if Q_p <= 0
%         Q_p = 0;
%         if P_rv < P_pa
%             Q_p_dot = 0;
%         end
%     end  
    
% Pulmonary Valve Flow
if P_pa >= P_rv
    Q_p = 0;
else
    Q_p = sqrt(P_rv - P_pa)/Rv;
end  
    
% Right atrial and ventricular volume
    % Systemic vein flow    
        Q_sv = (P_sv - P_ra) / RSV; %right side return, flow into RA
    % Right atrial volume ODE
        V_ra_dot = Q_sv - Q_t;     
    % Right ventricular volume ODE
        V_rv_dot = Q_t - Q_p;
     
% Pulmonary artery pressure and flow ODEs
    P_pa_dot = (Q_p - Q_pa)/CPA;
    Q_pa_dot = (P_pa - P_pt - Q_pa*RPA)/LPA;
    
% Pulmonary artery tree flow and pressure ODE
    Q_pt = (P_pt - P_pv) / RPT;
    P_pt_dot = (Q_pa - Q_pt) / CPT;    
   
% Pulmonary vein pressure ODE
    P_pv_dot = (Q_pt - Q_pv)/CPV;
    
% Systemic vein pressure ODE - had to place after Q_sv is defined
    P_sv_dot = (Q_st - Q_sv)/CSV;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdot = zeros(12,1);
xdot(1) = V_la_dot;    
xdot(2) = V_lv_dot;    
xdot(3) = P_sa_dot;  
xdot(4) = Q_sa_dot;
xdot(5) = P_st_dot; 
xdot(6) = P_sv_dot;    
xdot(7) = V_ra_dot;        
xdot(8) = V_rv_dot;        
xdot(9) = P_pa_dot;  
xdot(10)= Q_pa_dot;
xdot(11)= P_pt_dot; 
xdot(12)= P_pv_dot;
% xdot(13)= Q_m_dot;
% xdot(14)= Q_a_dot;
% xdot(15)= Q_t_dot;
% xdot(16)= Q_p_dot;


%% States
% xdot1 = [PSAdot;PSTdot;PSVdot;PPAdot;PPTdot;PPVdot;QSAdot;QPAdot;VLAdot;VLVdot;VRAdot;VRVdot];
% xdot2 = [QMIdot;QAOdot;QTValvedot;QPValvedot];
% xdot  = [xdot1;xdot2];

%% Outputs
% y(1) = PSA;
% y(2) = PLV;

% Stage of cardiac cycle
if Q_a ~= 0
    % Ejection
    stage = 1;
elseif Q_m ~= 0
    % Filling
    stage = 2;
elseif Q_a == 0 && Q_m == 0
    stage = 3;
end

y(1) = P_sa;
y(2) = P_lv;
y(3) = Q_a;
y(4) = Q_m;
y(5) = e_t_lv;
y(6) = Qvad;
y(7) = P_la;
y(8) = P_ra;
% Output R-wave trigger
if tn == 0
    rwt = 1;
else
    rwt = 0;
end
y(9) = rwt;
end














