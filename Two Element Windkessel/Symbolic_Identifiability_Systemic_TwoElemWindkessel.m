%==========================================================================
% Use Hermann & Krener nonlinear observability rank test, as well as
% nonlinear observability using the extended Lie derivatives approach
% (Villaverde) to investigate identifiability of systemic circulation 
% parameters of the two element Windkessel model.
%
% Available measurements - We will only use Pao and Qao for this part
%
% System states - Pao
% System inputs - Qao (Qvad + Qa)
%
% RESULTS (4-23-2020)
% We find that Pr and Cs are identifiable during all stages. We assume that
% we can measure Rsvr using MAP, CVP and mean CO.
%
%..........................................................................
% Suraj R Pawar, 4/23/2020
%==========================================================================

clear all; close all; clc;

%% Define symbols

    syms Pao;               % System states
    syms Cs Pr Rsvr;        % Unkonwn parameters (time invariant)   
    syms Qao;               % System input     
    
%% Symbolic nonlinear observability - Hermann and Krener approach

    % -------------  User input -------------------------------------------
        x = [Pao];                     % State vector
        theta = [Cs; Pr];          % Unknown parameters        
        y = [Pao];                     % Measurements
    %----------------------------------------------------------------------
    
    x_aug = [x; theta];
    num_params = length(theta);
    param_dynamics = zeros(num_params,1);
   
    %% Ejection
    
        f_ejection = [(1/Cs)*(Qao - (1/Rsvr)*(Pao + Pr));
                       param_dynamics];
        y_ejection = y;

        [r_ejection, ~] = func_NLO(x_aug,f_ejection, y_ejection)             
        
    %% Isovolumic 
        f_iso = [(1/Cs)*(Qao - (1/Rsvr)*(Pao + Pr));
                  param_dynamics];
        y_iso = y;
        
        [r_iso, ~] = func_NLO(x_aug,f_iso, y_iso);
            
    %% Filling
        f_filling = [(1/Cs)*(Qao - (1/Rsvr)*(Pao + Pr));
                      param_dynamics];
        y_filling = y;
        
        [r_filling, ~] = func_NLO(x_aug,f_filling, y_filling)
                
        
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
    syms dQao; % Consider upto the first order derivative of the input
    uvec = [Qao; dQao];
    
    %% Ejection
        r_ejection_villaverde = func_NLO_Villaverde(x_aug,f_ejection, y_ejection, uvec);
        
    %% Isovolumic
        r_isovolumic_villaverde = func_NLO_Villaverde(x_aug,f_iso, y_iso, uvec);
        
    %% Filling
        r_filling_villaverde = func_NLO_Villaverde(x_aug,f_filling, y_filling, uvec);
    
    %% Console Output
        fprintf('=== Nonlinear observability (Villaverde) ============\n\n');
        fprintf('Assuming that dQvad is non-zero\n');
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
    