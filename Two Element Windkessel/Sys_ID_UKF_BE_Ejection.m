%% B, Emax estimation using UKF and Two Element Windkessel model
%{
---------------------------- Description ----------------------------------
Uses Unscented Kalman Filter (UKF) to perform system identification of the
CVS using a two element model as the model basis. Assume that we know all
other parameters except for B, Emax. Estimation will be performed only
during isovolumic contraction and ejection phase. The process of determining 
the right stage of the cardiac cycle is as follows : 

1. Perform power spectrum analysis on aortic flow signal.
2. Choose a cutoff frequency that satisfactorily removes noise.
3. Apply lowpass filter to aortic flow.
4. Scan all samples, and start marking the samples when normalized time, i.e
   tn = mod(t,tc), is within 0 and dt. This marks the start of the R wave.
   Here tc is the time for one cardiac cycle. 
5. Keep scanning and marking the samples until aortic flow drops below some
   threshold value. Stop marking here. 
6. The samples you have marked correspond to isovolumic contraction and ejection.
7. Use these sample windows for UKF estimation.

------------------------- Assumptions -------------------------------------
1. We assume that we know Pr, A, SVR, Cao, Rv
2. UKF for B and Emax estimation was tested for various values of A
    - Seems to work within 10% accuracy when A is deviated +50% to -10%
3. SVR can be estimated using the formula
        SVR = (MAP - CVP) / CO
4. We could try to set a value of Rv from deltaP / CO
5. Pr and Cao may be estimable through the systemic circulation side
   estimation

-------------------------------- Verisons ---------------------------------
v1 : 5-26-2020, Suraj R Pawar
    - Initialize
v2 : 5-28-220, Suraj R Pawar
    - Implement versioning system to be able to extract parameters and set
    their dynamic equations based on the experiment being run
v3 : 6-9-2020, Suraj R Pawar
    - Added 'include_us.m' to add necessary folders and subfolders to
    MATLAB path
v4 : Suraj R Pawar, 6-11-2020
    - Tested after adding 'counts' as an extra input to
    func_handle_parameters
%}

close all; clear all; clc;
include_us;

%% User Inputs
    % Initial Conditions
    V0_guess = 5;                       % Guess for unstressed blood volume (mL)
    Vlv0 = 100;                         % Initial left ventricle volume (mL)
    Ps0 = 100;                          % Initial systemic pressure (mmHg)
    
    A_deviation =   0;                 % Percentage deviation from true A
    B_deviation =   0;                  % Percentage deviation from true B
    E_deviation =   0;                  % Percentage deviation from true E
    B0 = 0.01;
    Emax0 = 1;                          % Initial guess; Max LV Elastance
    
    p0 = diag([100, 50, 1e-3, 4]);      % Initial error covariance
        
    param_noise_std = 1*[1e-10; 1e-10]; % White noise standard deviation for parameters [B; Emax]
    
    % UKF parameters
    alpha = 1e-3;
    kappa = 0;
    beta = 2;                          % 2 is optimal for Gaussian distributions
    
    % What is the version of this experiment ? 
    experiment_versions;               % Contains description of versions
    version = 1;

%% Measurements and parameters
    data = noisy_twoelem_data;      % Load noisy two element data using this function
      
    % True value of parameters
    Cs_true = data.parameters(1);
    Rsvr_true = data.parameters(2);
    Pr_true = data.parameters(3);
    Ra_true = data.parameters(4);
    Rm_true = data.parameters(5);
    A_true = (1 + A_deviation/100)*data.parameters(6);    
    B_true = (1 + B_deviation/100)*data.parameters(7);
    E_true = (1 + E_deviation/100)*data.parameters(8);
    V0_true = data.parameters(9);
    HR_true = data.parameters(10);
    t_c_true = data.parameters(11);
    t_vc_true = data.parameters(12);    
    parameters = [Rsvr_true; Cs_true; Pr_true; Ra_true; Rm_true; HR_true; t_vc_true; t_c_true;...
                  A_true; B_true; E_true; V0_true]; % To be passed to UKF       

    % Filter aortic flow - will be used to detect ejection windows
    dt_original = data.dt;
    Fs = 1/dt_original;                         % Sampling frequency (Hz)
    Qa_original = data.Qa;
    Qa_filtered = lowpass(Qa_original, 30, Fs);
    
    % Interpolate all measurements to 1ms timing    
    dt = 0.001;
    tf = data.num_beats*data.parameters(11);   % Final time for simulation
    t_original = [0 : dt_original : tf];       % Time vector from measurement file
    t = [0: dt : tf];                          % Time vector with 1 ms time steps
    Plv_original = data.Plv;
    Ps_original = data.Ps;    
    Qvad_original = data.Qvad;
    Vbar_original = data.x(1,:) - V0_true;    
    
    Plv = interp1(t_original, Plv_original, t);
    Pao = interp1(t_original, Ps_original, t);
    Qa = interp1(t_original, Qa_original, t);
    Qa_filtered = interp1(t_original, Qa_filtered, t);
    Qvad = interp1(t_original, Qvad_original, t);
    Vbar = interp1(t_original, Vbar_original, t);   % Vbar = Vlv - V0    
    y = [Plv; Pao; Qa];    
    
    % Noise and parameters for UKF
    %{
    Further comments : 
    This is the white Gaussian noise strength for the continuous time
    model. We need to convert this into Ito's form and integrate using
    Euler integration. When we do convert it to Ito's form inside the
    function, the variance of the Brownian motion is q*dt
    %}
    q = diag([data.process_noise_std.^2; param_noise_std.^2]);  % Process noise covariance
    r = diag(data.meas_noise_std.^2);                           % Measurement noise covariance
    ukf_params = [alpha; kappa; beta];
    x0 = [Vlv0 - V0_guess; Ps0];                                % First state is Vbar = Vlv - V0
    theta0 = [B0; Emax0];
    
%% UKF 
    waitflag = 1;
    [xhat, yhat, Paug, teject] = func_TwoElem_SysID_UKF_BE_Ejection(t, y, Qvad, x0, theta0, p0,...
                                                                   q, r, parameters, ukf_params,...
                                                                   version, Qa_filtered, waitflag);

%% Figures
    figure;   
    
    subplot(3,2,1); % Vbar
    hold on;
    plot(t,Vbar);
    plot(teject,xhat(1,:),'r.');
    hold off;
    legend({'Measured','Estimated'},'Orientation','horizontal');
    title('Vbar (mL)');
    xlabel('Time (s)');
    
    subplot(3,2,2); % Plv
    hold on;
    plot(t,Plv);
    plot(teject,yhat(1,:),'r.');
    hold off;
    legend({'Measured','Estimated'},'Orientation','horizontal');
    title('Plv (mmHg)');
    xlabel('Time (s)');
    
    subplot(3,2,3); % Ps
    hold on;
    plot(t,Pao);
    plot(teject,xhat(2,:),'r.');
    hold off;
    legend({'Measured','Estimated'},'Orientation','horizontal');
    title('Pao (mmHg)');
    xlabel('Time (s)');    
    
    subplot(3,2,4); % B
    hold on;
    plot(t,B_true*ones(size(t)),'--k');
    plot(teject,xhat(3,:));
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('B (1/mL)');
    xlabel('Time (s)');
    
    subplot(3,2,5); % Emax
    hold on;
    plot(t,E_true*ones(size(t)),'--k');
    plot(teject,xhat(4,:));
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('Emax (mmHg/mL)');
    xlabel('Time (s)');