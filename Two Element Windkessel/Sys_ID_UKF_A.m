%% A estimation using UKF and Two Element Windkessel model
%{
---------------------------- Description ----------------------------------
Uses Unscented Kalman Filter (UKF) to perform system identification of the
CVS using a two element model as the model basis. Assume that we know all
other parameters except for A. Estimation will be performed only
during iso expansion and filling phase. The process of determining 
the right stage of the cardiac cycle is as follows : 

1. Perform power spectrum analysis on aortic flow signal.
2. Choose a cutoff frequency that satisfactorily removes noise.
3. Apply lowpass filter to aortic flow.
4. Scan all samples, and start marking the samples when normalized time, i.e
   tn = mod(t,tc) is within 0 and dt. This marks the start of the R wave.
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
v3 : Suraj R pawar, 6-9-2020
    - Add 'include_us' file
%}

close all; clear all; clc;
include_us;

%% User Inputs
    % Initial Conditions
    V0_guess = 5;                       % Guess for unstressed blood volume
    Vlv0 = 10;                         % Initial left ventricle volume (mL)
    Ps0 = 100;                          % Initial systemic pressure (mmHg)
    
    A_deviation =   0;                  % Percentage deviation from true A
    B_deviation =   0;                  % Percentage deviation from true B
    E_deviation =   0;                  % Percentage deviation from true E
    A0 = 0.01;    
    
    p0 = diag([100, 25, 1e-1]);      % Initial error covariance
        
    param_noise_std = 1*[1e-15]; % White noise standard deviation for parameters [A]
    
    % UKF parameters
    alpha = 1e-3;
    kappa = 0;
    beta = 2;                          % 2 is optimal for Gaussian distributions
    
    % What is the version of this experiment ? 
    experiment_versions;               % Contains description of versions
    version = 2;
        
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

    % Filter Plv - will be used to detect ejection windows
    dt_original = data.dt;
    Fs = 1/dt_original;                         % Sampling frequency (Hz)
    Plv_original = data.Plv;
    Plv_filtered = lowpass(Plv_original, 4, Fs);
    
    % Interpolate all measurements to 1ms timing    
    dt = 0.001;
    tf = data.num_beats*data.parameters(11);   % Final time for simulation
    t_original = [0 : dt_original : tf];       % Time vector from measurement file
    t = [0: dt : tf];                          % Time vector with 1 ms time steps    
    Ps_original = data.Ps;    
    Qvad_original = data.Qvad;
    Qa_original = data.Qa;
    Vbar_original = data.x(1,:) - V0_true;    
    
    Plv = interp1(t_original, Plv_original, t);
    Plv_filtered = interp1(t_original, Plv_filtered, t);
    Pao = interp1(t_original, Ps_original, t);
    Qa = interp1(t_original, Qa_original, t);
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
    theta0 = [A0];
    
%% UKF 
    [xhat, yhat, Paug, tselect] = func_TwoElem_SysID_UKF_A(t, y, Qvad, x0, theta0, p0,...
                                                                   q, r, parameters, ukf_params,...
                                                                   version, Plv_filtered);

%% Figures
    figure;   
    
    subplot(2,2,1); % Vbar
    hold on;
    plot(t,Vbar);
    plot(tselect,xhat(1,:),'r.');
    hold off;
    legend({'Measured','Estimated'},'Orientation','horizontal');
    title('Vbar (mL)');
    xlabel('Time (s)');
    
    subplot(2,2,2); % Plv
    hold on;
    plot(t,Plv);
    plot(tselect,yhat(1,:),'r.');
    hold off;
    legend({'Measured','Estimated'},'Orientation','horizontal');
    title('Plv (mmHg)');
    xlabel('Time (s)');
    
    subplot(2,2,3); % Ps
    hold on;
    plot(t,Pao);
    plot(tselect,xhat(2,:),'r.');
    hold off;
    legend({'Measured','Estimated'},'Orientation','horizontal');
    title('Pao (mmHg)');
    xlabel('Time (s)');    
    
    subplot(2,2,4); % A
    hold on;
    plot(t,A_true*ones(size(t)),'--k');
    plot(tselect,xhat(3,:));
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('A (mmHg)');
    xlabel('Time (s)');        