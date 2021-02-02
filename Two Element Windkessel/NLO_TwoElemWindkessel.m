%==========================================================================
% Use Hermann & Krener nonlinear observability rank test, as well as
% nonlinear observability using the extended Lie derivatives approach
% (Villaverde) to investigate identifiability of parameters of the two
% element Windkessel model. 
%
% Assumptions - 
% 1. Plv equation as given in Gohean 2013
% 2. (Vlv - V0) taken as a new state, and this avoids counting V0 as a
%    separate parameter
% 3. Assume that SVR can be measured using (mean PAO - mean CVP)/mean CO.
%    We already have measurment signals required in that equation
% 4. For NLO formulation, we are not allowed to keep time-varying
%    parameters. Since we known the value of en, we can treat it as a
%    time-varying input and calculate observability rank as per Villaverde.
%
% Available measurements - Plv, Pao, Qao (aortic flow, does not include VAD
% flow), Qvad (Qvad to be treated as input to the system).
%
% Possible measurements that could be included in future experiments - Vlv
%
% System states - Vlv, Pao
% 
%..........................................................................
% Suraj R Pawar, 4/17/2020
%==========================================================================

clear all; close all; clc;

%% Define symbols

    syms Vlv Pao;               % System states
    syms Rv Cao Pr;    % Unkonwn parameters (time invariant), assume Ra = Rm = Rv
    syms Rsvr;                  % Known parameters
    syms Qvad Plv;               % System input     
    
%% Symbolic nonlinear observability - Hermann and Krener approach

    % -------------  User input -------------------------------------------
        x = [Vlv; Pao];                                 % State vector
        theta = [Cao; Rv; Pr];              % Unknown parameters                
        Qa = (Plv - Pao)/Rv;                            % Equation for aortic flow
        y = [Pao; Qa];                             % Measurements. Qa will be removed for iso and filling stages       
    %----------------------------------------------------------------------
    
    x_aug = [x; theta];
    num_params = length(theta);
    param_dynamics = zeros(num_params,1);
   
    %% Ejection
    
        f_ejection = [-Qvad - (Plv - Pao)/Rv;
                       (1/Cao)*(Qvad - (1/Rsvr)*(Pao) + (Plv - Pao)/Rv);
                       param_dynamics];
        y_ejection = y;
%         y_ejection = Plv - Pao;
        [r_ejection, ~] = func_NLO(x_aug,f_ejection, y_ejection);        
        
    %% Isovolumic 
        f_iso = [-Qvad;
                  (1/Cao)*(Qvad - (1/Rsvr)*(Pao));
                  param_dynamics];
        y_iso = [y(1); 0];
%         y_iso = Plv - Pao;
        
        [r_iso, ~] = func_NLO(x_aug,f_iso, y_iso);
            
    %% Filling
        f_filling = [-Qvad + (Pr - Plv)/Rv;
                      (1/Cao)*(Qvad - (1/Rsvr)*(Pao));
                      param_dynamics];
        y_filling = [y(1); 0];
%         y_filling = Plv - Pao;
        [r_filling, ~] = func_NLO(x_aug,f_filling, y_filling);           
        
    %% Console Output
        fprintf('=== Nonlinear observability (Hermann and Krener) =====\n\n');
        fprintf('EJECTION \n');
        fprintf('Number of states in augmented state vector: %d \n', length(x_aug));
        fprintf('Nonlinear observability rank: %d \n', r_ejection);
        fprintf('\n\n');
        
        fprintf('ISOVOLUMIC \n');
        fprintf('Number of states in augmented state vector: %d \n', length(x_aug));
        fprintf('Nonlinear observability rank: %d \n', r_iso);
        fprintf('\n\n');
        
        fprintf('FILLING \n');
        fprintf('Number of states in augmented state vector: %d \n', length(x_aug));
        fprintf('Nonlinear observability rank: %d \n', r_filling);
        fprintf('\n\n');
        
%% Symbolic nonlinear observability - Villaverde approach
    syms dQvad dPlv; % Consider upto the first order derivative of the input
    uvec = [Qvad, dQvad; Plv, dPlv];
    
    %% Ejection
        r_ejection_villaverde = func_NLO_Villaverde(x_aug,f_ejection, y_ejection, uvec);
        
    %% Isovolumic
        r_isovolumic_villaverde = func_NLO_Villaverde(x_aug,f_iso, y_iso, uvec);
        
    %% Filling
        r_filling_villaverde = func_NLO_Villaverde(x_aug,f_filling, y_filling, uvec);
    
    %% Console Output
        fprintf('=== Nonlinear observability (Villaverde) ============\n\n');
        fprintf('Assuming that dQvad and dPlv are non-zero\n');
        fprintf('EJECTION \n');
        fprintf('Number of states in augmented state vector: %d \n', length(x_aug));
        fprintf('Nonlinear observability rank: %d \n', r_ejection_villaverde);
        fprintf('\n\n');
        
        fprintf('ISOVOLUMIC \n');
        fprintf('Number of states in augmented state vector: %d \n', length(x_aug));
        fprintf('Nonlinear observability rank: %d \n', r_isovolumic_villaverde);
        fprintf('\n\n');
        
        fprintf('FILLING \n');
        fprintf('Number of states in augmented state vector: %d \n', length(x_aug));
        fprintf('Nonlinear observability rank: %d \n', r_filling_villaverde);
        fprintf('\n\n');
    