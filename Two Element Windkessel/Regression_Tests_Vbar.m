% Regression tests for system identification
%{
----------------------- Description ---------------------------------------
Use this script to run regression test where we vary the initial Vbar0. 
The UKF is designed based on the Two Element Windkessel model with a square 
root law for valve flows. 

-------------------------- Versions ---------------------------------------
%{
v1 : Suraj R Pawar, 7-24-2020
    - Initialize
%}
%}

clear all; close all; clc;
include_us;

%% ------------------------- User Inputs ----------------------------------
    
    % Experiment
    Vset = 225;
    delta = 10;
    Vbar0_lower = (1 - delta/100)*Vset;
    Vbar0_upper = (1 + delta/100)*Vset;
    Vbar0_increment = 3;
    
    num_cycles = 30;            % Number of heart cycles
    vad = 1;
    
    % Logging
    filename = 'Results/Regression_Tests.xlsx';
    logging = 1;        
    
    parameters = [0.18, 0.016, 0.3, 0.65, 1.075, 0.0025, 90];       % A, B, E, Cs, Rsvr, Rv, HR
    heart_condition = 2;                                            % 1 = Healthy, 2 = Heart Failure
    
    override_tvc = 0;                                               % 0 = Calculate tvc, 1 = Manually set tvc 
    tvc_manual = 0.6;   
    
    % Initial Guesses
    initial_guesses.A0 = 0.01;
    initial_guesses.B0 = 0.01;
    initial_guesses.Emax0 = 1;
    %initial_guesses.Vbar0 = 150;        % 150 for healthy, 250 for heart failure
    
    % Filtering for Qa signal
    Qa_filter.lowpass = 30;
    Qa_filter.upper = 50;
    Qa_filter.lower = 0;
    
    % Known / Approximated parameters
    Pr_true = 14;    % (mmHg), 3 for healthy, 14 for heart failure
    V0_true = 5;    % (mL)
        
    % Noise statistics   
    pressure_noise_process = 5;  
    flow_noise_process = 1;   % Not used in UKF weightings
    volume_noise_process = 0.5;
    
    pressure_noise_meas = 0.8; 
    flow_noise_meas = 0.5;
    volume_noise_meas = 0.8;    % Not used in UKF weightings
    
    % Plotting parameters
    decimation = 10;
    figurestyle = 2;            % 1 = separate figures, 2 = all together
    
%% ----------------------------- Prepare parameters -----------------------
    Vbar0 = [Vbar0_lower : Vbar0_increment : Vbar0_upper];
    num_trials = numel(Vbar0);
    
    process_noise = [pressure_noise_process; volume_noise_process; flow_noise_process];
    meas_noise = [pressure_noise_meas; flow_noise_meas; volume_noise_meas];
    
%% Regression test
    
    % Memory Allocation    
    parameters_est = zeros(num_trials,6);
    RMSE = zeros(4,num_trials);
    
    for i = 1:num_trials
        fprintf('Trial number %d of %d\n',i,num_trials);                
        
        % Set initial Vbar0
        initial_guesses.Vbar0 = Vbar0(i);
        
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
        
        % RMSEs
        RMSE(:,i) = Variation_Compare(parameters_est(i,:), parameters, signals, noisy_signals, Vbar0(i), Pr_true, num_cycles);
    end
    
%% Data Logging
    if logging == 1        
        A_est = parameters_est(:,1);        
        B_est = parameters_est(:,2);        
        E_est = parameters_est(:,3);
        Cs_est = parameters_est(:,4);        
        svr_est = parameters_est(:,5);
        Rv_est = parameters_est(:,6);
                
        T = [Vbar0.', A_est, B_est, E_est, Cs_est, svr_est, Rv_est, RMSE.'];
              
        xlswrite(filename,T);
        fprintf('Data written to file .. \n');
    end
    
%% Figures
    
    pa = polyfit(Vbar0.',A_est,2);
    pb = polyfit(Vbar0.',B_est,2);
    pe = polyfit(Vbar0.',E_est,1);
    pc = polyfit(Vbar0.',Cs_est,1);
    pr = polyfit(Vbar0.',Rv_est,1);
    
    Vbar_xrange = [0.9*Vbar0_lower : Vbar0_increment : 1.1*Vbar0_upper];
    ya = polyval(pa,Vbar_xrange);
    yb = polyval(pb,Vbar_xrange);
    ye = polyval(pe,Vbar_xrange);
    yc = polyval(pc,Vbar_xrange);
    yr = polyval(pr,Vbar_xrange);
    
    figure;
    hold on;
    plot(Vbar0,A_est,'ok');
    plot(Vbar_xrange,ya,'-','Color',[0.6,0.6,0.6]);         
    plot(Vbar0,B_est,'xk'); hold on;
    plot(Vbar_xrange,yb,'-','Color',[0.6,0.6,0.6]);
    hold off;    
    %xlim([80 220]);
    %ylim([-0.01 0.08]);