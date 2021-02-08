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
v2 : Suraj R Pawar, 2-7-2021
    - Animal experiment tests with new two element model
%}
v3 : Suraj R Pawar, 2-8-2021
    - Pr now taken as end diastolic PLV. This is scanned from PLV waveform
    when computational model generates this signal
%}

clear all; close all; clc;
include_us;

%% ------------------------- User Inputs ----------------------------------
    
    % Experiment
    num_trials = 11;
    num_cycles = 30;
    vad = 1;
    
    % Logging
    filename = 'Results/Regression_Tests.xlsx';
    logging = 1;
    
    % Upper and lower Limits for parameters
    uA = 0.216; 
    lA = 0.144;
    
    uB = 0.0192;
    lB = 0.0128;
    
    uE = 0.36;
    lE = 0.24;
    
    uCs = 0.78;
    lCs = 0.52;
    
    uRsvr = 1.1;
    lRsvr = 1;
    
    uRv = 0.003;
    lRv = 0.002;
    
    uHR = 90;
    lHR = 90;
    
    override_parameters = 0;                                        % If you want to manually set the parameters
    parameters_set = [0.18, 0.016, 0.3, 0.65, 1.075, 0.0025, 90];   % A, B, E, Cs, Rsvr, Rv, HR
    %parameters_set = [0.0312    0.0572    3.3819    1.2754    0.9578    0.0025   80.0000];   % A, B, E, Cs, Rsvr, Rv, HR
    heart_condition = 2;                                            % 1 = Healthy, 2 = Heart Failure
    
    override_tvc = 0;                                               % 0 = Calculate tvc, 1 = Manually set tvc 
    tvc_manual = 0.6; %0.6   
    
    % Initial Guesses
    initial_guesses.A0 = 0.1;
    initial_guesses.B0 = 0.01;
    initial_guesses.Emax0 = 0.1;
    initial_guesses.Vbar0 = 250;        % 150 for healthy, 250 for heart failure
    
    % Filtering for Qa signal
    Qa_filter.lowpass = 30;
    Qa_filter.upper = 50;
    Qa_filter.lower = 0;
    
    % Known / Approximated parameters
    %Pr_true = 3;    % (mmHg), 3 for healthy, 14 for heart failure
    Pr_true = 3;    % (mmHg), 3 for healthy, 14 for heart failure
    V0_true = 5;    % (mL)
        
    % Noise statistics   
    pressure_noise_process = 5;  
    flow_noise_process = 1;   % Not used in UKF weightings
    volume_noise_process = 0.5; %0.5
    
    pressure_noise_meas = 0.8; 
    flow_noise_meas = 0.5;
    volume_noise_meas = 0.8;    % Not used in UKF weightings
    
    % Plotting parameters
    decimation = 10;
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
        [signals, Pr_true] = comp_data_generator(parameters, num_cycles, vad, heart_condition);
        
        % Load data from animal experiment
        %{
        data = noisy_twoelem_data('sync_1_1_SysID.mat');    
        num_elems = numel(data.Plv);
        zv = zeros(num_elems,1);
        signals = [data.Plv.', data.Ps.', zv, zv, zv, data.Qa.', data.Qvad.',...
                   data.t.'];
        
        parameters = [data.parameters(6),data.parameters(7),data.parameters(8),...
                      data.parameters(1), data.parameters(2), data.parameters(4),...
                      data.parameters(10)];
        %}
        
        % Add measurement noise
        noisy_signals = add_meas_noise(signals,meas_noise);
        %noisy_signals = signals;
        
        % Estimate SVR
        svr = svr_estimator(noisy_signals);
        %svr = 0.7684;   % (mmHg.s/mL) This is manually measured from animal experiment dataset
        
        % Estimate Cs and Rv
        if override_tvc == 0
            t_vc_true = (550 - 1.75*parameters(7))/1000;
        elseif override_tvc == 1
            t_vc_true = tvc_manual;
        end
        
        [Cs, Rv, y1] = CsRv_estimator(noisy_signals, parameters, Pr_true, V0_true, t_vc_true, process_noise, meas_noise);
        if Cs < 0
            fprintf('Bad Run... \n');
            %return;
        end
        
        % Estimate A, B, Emax                
        initial_guesses.Ps0 = noisy_signals(1,2);
        %initial_guesses.Vbar0 = noisy_signals(1,5);
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
                
        T = [A_true,A_est,B_true, B_est, E_true, E_est, ...
                  Cs_true, Cs_est, svr_true, svr_est, Rv_true, Rv_est];
              
        xlswrite(filename,T);
        fprintf('Data written to file .. \n');
    end