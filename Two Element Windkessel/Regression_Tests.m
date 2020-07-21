% Regression tests for system identification
%{
----------------------- Description ---------------------------------------
Use this script to run all sorts of regression tests with the system
identification code based on UKF. The UKF is designed based on the Two
Element Windkessel model with a square root law for valve flows. 

-------------------------- Versions ---------------------------------------
%{
v1 : Suraj R Pawar, 7-20-2020
    - Initialize
%}
%}

clear all; close all; clc;
include_us;

%% ------------------------- User Inputs ----------------------------------
    
    % Experiment
    num_trials = 1;
    num_cycles = 30;
    vad = 1;
    
    % Logging
    filename = 'Results/Regression_Tests.xlsx';
    logging = 0;
    
    % Upper and lower Limits for parameters
    uA = 0.2; 
    lA = 0.15;
    
    uB = 0.018;
    lB = 0.012;
    
    uE = 0.35;
    lE = 0.28;
    
    uCs = 0.7;
    lCs = 0.6;
    
    uRsvr = 1.08;
    lRsvr = 1.05;
    
    uRv = 0.003;
    lRv = 0.002;
    
    uHR = 82;
    lHR = 80;
    
    override_parameters = 1;                                        % If you want to manually set the parameters
    parameters_set = [0.18, 0.016, 0.3, 0.65, 1.075, 0.0025, 90];   % A, B, E, Cs, Rsvr, Rv, HR
    heart_condition = 2;                                            % 1 = Healthy, 2 = Heart Failure
    
    override_tvc = 0;                                               % 0 = Calculate tvc, 1 = Manually set tvc 
    tvc_manual = 0.6;   
    
    % Initial Guesses
    initial_guesses.A0 = 0.01;
    initial_guesses.B0 = 0.01;
    initial_guesses.Emax0 = 1;
    initial_guesses.Vbar0 = 280;        % 150 for healthy, 250 for heart failure
    
    % Filtering for Qa signal
    Qa_filter.lowpass = 15;
    Qa_filter.upper = 50;
    Qa_filter.lower = 0;
    
    % Known / Approximated parameters
    Pr_true = 14;    % (mmHg), 3 for healthy, 14 for heart failure
    V0_true = 5;    % (mL)
        
    % Noise statistics   
    pressure_noise_process = 5;  
    flow_noise_process = 0.8;   % Not used in UKF weightings
    volume_noise_process = 1;
    
    pressure_noise_meas = 0.8; 
    flow_noise_meas = 0.1;
    volume_noise_meas = 0.8;    % Not used in UKF weightings
    
    % Plotting parameters
    decimation = 5;
    figurestyle = 2;            % 1 = separate figures, 2 = all together
    
%% ----------------------------- Prepare parameters -----------------------
    upper_limits = [uA; uB; uE; uCs; uRsvr; uRv; uHR];
    lower_limits = [lA; lB; lE; lCs; lRsvr; lRv; lHR];
    
    process_noise = [pressure_noise_process; volume_noise_process; flow_noise_process];
    meas_noise = [pressure_noise_meas; flow_noise_meas; volume_noise_meas];
    
%% Regression test
    
    % Memory Allocation
    parameters_true = zeros(num_trials,7); % Also includes heart rate
    parameters_est = zeros(num_trials,6);
    
    for i = 1:num_trials
        fprintf('Trial number %d of %d\n',i,num_trials);
        
        % Parameter Generator
        if override_parameters == 1
            parameters_true(i,:) = parameters_set;
        else
            parameters_true(i,:) = param_generator(upper_limits, lower_limits);
        end
        
        parameters = parameters_true(i,:);
        
        % Generate Measurements
        signals = comp_data_generator(parameters, num_cycles, vad, heart_condition);

        % Add measurement noise
        noisy_signals = add_meas_noise(signals,meas_noise);

        % Estimate SVR
        svr = svr_estimator(noisy_signals);

        % Estimate Cs and Rv
        if override_tvc == 0
            t_vc_true = (550 - 1.75*parameters(7))/1000;
        elseif override_tvc == 1
            t_vc_true = tvc_manual;
        end
        
        [Cs, Rv, y1] = CsRv_estimator(noisy_signals, parameters, Pr_true, V0_true, t_vc_true, process_noise, meas_noise);
        if Cs < 0
            fprintf('Aborting... \n');
            return;
        end
        
        % Estimate A, B, Emax                
        initial_guesses.Ps0 = noisy_signals(1,2);
        [A, B, Emax, y2] = ABE_estimator(noisy_signals, parameters, Pr_true, Cs, Rv, svr, t_vc_true, process_noise, meas_noise, Qa_filter, initial_guesses, num_cycles);
        
        % Collect estimated parameters
        parameters_est(i,:) = [A, B, Emax, Cs, svr, Rv];
        
        % Figures
        if num_trials == 1                        
            figures(signals, noisy_signals, y1, y2, V0_true, decimation, figurestyle, parameters_est(i,:));
        end
    end
    
%% Data Logging
    if logging == 1
        A_true = parameters_true(:,1);
        A_est = parameters_est(:,1);
        
        B_true = parameters_true(:,2);
        B_est = parameters_est(:,2);
        
        E_true = parameters_true(:,3);
        E_est = parameters_est(:,3);
        
        Cs_true = parameters_true(:,4);
        Cs_est = parameters_est(:,4);
        
        svr_true = parameters_true(:,5);
        svr_est = parameters_est(:,5);
        
        Rv_true = parameters_true(:,6);
        Rv_est = parameters_est(:,6);
        
        T = table(A_true,A_est,B_true, B_est, E_true, E_est, ...
                  Cs_true, Cs_est, svr_true, svr_est, Rv_true, Rv_est);
              
        writetable(T,filename);
    end