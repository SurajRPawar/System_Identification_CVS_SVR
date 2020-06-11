%% Cs and Pr estimation using UKF and Two Element Windkessel model
%{
---------------------------- Description ----------------------------------
Uses Unscented Kalman Filter (UKF) to perform system identification of the
CVS using a two element model as the model basis. Assume that we know all
other parameters except for Cs and Pr. Estimation of these parameters can 
be performed during all stages of the cardiac cycle. 

------------------------- Assumptions -------------------------------------
1. We assume that we know Rsvr, A, B, Emax, Rv
2. UKF can be used to estimation A, B, Emax
3. Rv can be approximated as the mean dP / mean Qa
4. SVR can be estimated using the formula
        SVR = (MAP - CVP) / CO

-------------------------------- Verisons ---------------------------------
v1 : 6-11-2020, Suraj R Pawar
    - Initialize
%}

close all; clear all; clc;
include_us;

%% User Inputs
    % Initial Conditions    
    Ps0 = 120;                      % Initial systemic pressure (mmHg)        
    
    Cs0 = 5;                        % Initial guess for systemic compliance (mL/mmHg)
    Pr0 = 10;                       % Initial guess for Pr (mmHg)
    
    p0 = diag([250, 25, 1000]);    % Initial error covariance
        
    param_noise_std = 1*[1e-8; 1e-8];    % White noise standard deviation for parameters [A]
    
    % UKF parameters
    alpha = 1e-3;
    kappa = 0;
    beta = 2;                       % 2 is optimal for Gaussian distributions
    
    % What is the version of this experiment ? 
    experiment_versions;            % Contains description of versions
    version = 4;
        
%% Measurements and parameters
    data = noisy_twoelem_data;      % Load noisy two element data using this function
      
    % True value of parameters
    Cs_true = data.parameters(1);
    Rsvr_true = data.parameters(2);
    Pr_true = data.parameters(3);
    Ra_true = data.parameters(4);
    Rm_true = data.parameters(5);
    A_true = data.parameters(6);    
    B_true = data.parameters(7);
    E_true = data.parameters(8);
    V0_true = data.parameters(9);
    HR_true = data.parameters(10);
    t_c_true = data.parameters(11);
    t_vc_true = data.parameters(12);    
    parameters = [Rsvr_true; Cs_true; Pr_true; Ra_true; Rm_true; HR_true; t_vc_true; t_c_true;...
                  A_true; B_true; E_true; V0_true]; % To be passed to UKF       
    
    
    % Interpolate all measurements to 1ms timing    
    dt = 0.001;
    dt_original = data.dt;
    tf = data.num_beats*data.parameters(11);   % Final time for simulation
    t_original = [0 : dt_original : tf];       % Time vector from measurement file
    t = [0: dt : tf];                          % Time vector with 1 ms time steps    
    Ps_original = data.Ps;    
    Qvad_original = data.Qvad;
    Qa_original = data.Qa;    
         
    Pao = interp1(t_original, Ps_original, t);
    Qa = interp1(t_original, Qa_original, t);    
    Qvad = interp1(t_original, Qvad_original, t);    
    Qao = Qa + Qvad;
    
    y = [Pao];    
    
    % Noise and parameters for UKF
    %{
    Further comments : 
    This is the white Gaussian noise strength for the continuous time
    model. We need to convert this into Ito's form and integrate using
    Euler integration. When we do convert it to Ito's form inside the
    function, the variance of the Brownian motion is q*dt
    %}
    Ps_noise_std = data.process_noise_std(2);
    meas_noise_std = data.meas_noise_std(2);
    q = diag([Ps_noise_std.^2; param_noise_std.^2]);  % Process noise covariance
    r = diag(meas_noise_std.^2);                      % Measurement noise covariance
    ukf_params = [alpha; kappa; beta];
    x0 = [Ps0];                                       
    theta0 = [Cs0; Pr0];
    
%% UKF 
    fprintf('Beginning UKF estimation \n'); tic;
    waitflag = 1;
    [xhat, yhat, Paug] = func_TwoElem_SysID_UKF_CsPr(t, y, Qao, x0, theta0, p0,...
                                                                   q, r, parameters, ukf_params,...
                                                                   version, waitflag);
    fprintf('UKF estimation finished in %.2f seconds\n', toc);

%% Figures
    figure;   
    Psplot = subplot(3,1,1);
    Csplot = subplot(3,1,2);
    Prplot = subplot(3,1,3);
    ax = [];
    
    axes(Psplot); % Ps
    hold on;
    plot(t,Pao,'Color',[0.6,0.6,0.6],'LineWidth',2);
    plot(t,xhat(1,:),'b');
    hold off;
    title('Pao (mmHg)');
    legend({'Measured','Estimated'},'Orientation','horizontal');        
    ax = [ax, gca];
            
    axes(Csplot); % Cs
    hold on;
    plot(t,Cs_true*ones(size(t)),'--k');
    plot(t,xhat(2,:));
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('Cs (mL/mmHg)');         
    ax = [ax, gca];
    
    axes(Prplot); % Pr
    hold on;
    plot(t,Pr_true*ones(size(t)),'--k');
    plot(t,xhat(3,:));
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('Pr (mmHg)');
    xlabel('Time (s)');        
    ax = [ax, gca];
    
    linkaxes(ax,'x');