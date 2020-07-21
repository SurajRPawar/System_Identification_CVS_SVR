%% Cs and Rv estimation using UKF and Two Element Windkessel model
%{
---------------------------- Description ----------------------------------
Uses Unscented Kalman Filter (UKF) to perform system identification of the
CVS using a two element model as the model basis. Assume that we know all
other parameters except for Cs and Rv. Estimation of these parameters can 
be performed during all stages of the cardiac cycle. 

------------------------- Assumptions -------------------------------------
1. We assume that we know Rsvr, A, B, Emax, Pr
2. UKF can be used to estimation A, B, Emax
3. SVR can be estimated using the formula
        SVR = (MAP - CVP) / CO
4. Pr can be approximated from the Plv value during filling stage

-------------------------------- Verisons ---------------------------------
v1 : 6-22-2020, Suraj R Pawar
    - Initialize
%}

close all; clear all; clc;
include_us;

%% ---------------------- User Inputs -------------------------------------
    % Data file to load
    datafilename = 'noisy_data_comp_hf.mat';
    
    % Set known values
    Pr_true = 3.7;                        % (mmHg)
    
    % Initial Conditions    
    Ps0 = 65;                          % Initial systemic pressure (mmHg)            
    Cs0 = 1;                            % Initial guess for systemic compliance (mL/mmHg)    
    Rv0 = 0.001;    
    
    % Initial covariances
    p0 = diag([5, 1, 1e-5]);           % Initial error covariance
        
    % Process noise terms
    data.process_noise_std(2) = 1;      % If you want to override default value
    param_noise_std = 1*[1e-10; 1e-10]; % White noise standard deviation for parameters [A]
    
    % UKF parameters
    alpha = 1e-3;
    kappa = 0;
    beta = 2;                           % 2 is optimal for Gaussian distributions
    
    % What is the version of this experiment ? 
    experiment_versions;                % Contains description of versions
    version = 11;
        
    % GUI progress bar
    waitflag = 0;
    
%% ---------------------- Measurements and parameters ---------------------
    data = noisy_twoelem_data(datafilename);      % Load noisy two element data using this function
      
    % True value of parameters
    Cs_true = data.parameters(1);
    Rsvr_true = data.parameters(2);    
    Ra_true = data.parameters(4);
    Rm_true = data.parameters(5);
    A_true = data.parameters(6);    
    B_true = data.parameters(7);
    E_true = data.parameters(8);
    V0_true = data.parameters(9);
    HR_true = data.parameters(10);
    t_c_true = data.parameters(11);
    t_vc_true = 0.41;%data.parameters(12);    
    parameters = [Rsvr_true; Cs_true; Pr_true; Ra_true; Rm_true; HR_true; t_vc_true; t_c_true;...
                  A_true; B_true; E_true; V0_true]; % To be passed to UKF       
    
    
    % Interpolate all measurements to 1ms timing    
    dt = 0.001;
    dt_original = data.dt;
    tf = data.num_beats*data.parameters(11);   % Final time for simulation
    t_original = [0 : dt_original : tf];       % Time vector from measurement file
    t = [0: dt : tf];                          % Time vector with 1 ms time steps    
    Plv_original = data.Plv;
    Ps_original = data.Ps;    
    Qvad_original = data.Qvad;
    Qa_original = data.Qa;    
         
    Plv = interp1(t_original, Plv_original, t, 'linear', 'extrap');
    Pao = interp1(t_original, Ps_original, t, 'linear', 'extrap');
    Qa = interp1(t_original, Qa_original, t, 'linear', 'extrap');
    Qvad = interp1(t_original, Qvad_original, t, 'linear', 'extrap');   
    
    y = [Pao; Qa];
    u = [Plv; Qvad];
    
    % Noise and parameters for UKF
    %{
    Further comments : 
    This is the white Gaussian noise strength for the continuous time
    model. We need to convert this into Ito's form and integrate using
    Euler integration. When we do convert it to Ito's form inside the
    function, the variance of the Brownian motion is q*dt
    %}        
    meas_noise_std = data.meas_noise_std(2:3);
    Ps_noise_std = data.process_noise_std(2);
    q = diag([Ps_noise_std.^2; param_noise_std.^2]);  % Process noise covariance
    r = diag(meas_noise_std.^2);                      % Measurement noise covariance
    ukf_params = [alpha; kappa; beta];
    x0 = [Ps0];                                       
    theta0 = [Cs0; Rv0];
    
%% ----------------------------- UKF --------------------------------------
    fprintf('Beginning UKF estimation \n'); tic;    
    [xhat, yhat, Paug] = func_TwoElem_SysID_UKF_CsRv(t, y, u, x0, theta0, p0,...
                                                                   q, r, parameters, ukf_params,...
                                                                   version, waitflag);
    fprintf('UKF estimation finished in %.2f seconds\n', toc);

%% -------------------------- Console output ------------------------------
    Csmean = mean(xhat(2,[end-50:end]));   
    Rvmean = mean(xhat(3,[end-50:end]));   
    fprintf('Estimated Cs : %.3f mL/mmHg\n', Csmean);    
    fprintf('Estimated Rv : %.3f mmHg/mL/s\n', Rvmean);   
    
%% ------------------------------ Figures ---------------------------------
    figure;   
    numrows = 3;
    numcols = 1;
    Psplot = subplot(numrows,numcols,1);
    Csplot = subplot(numrows,numcols,2);
    Rvplot = subplot(numrows,numcols,3);
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
    
    axes(Rvplot); % Rv
    hold on;
    plot(t,Ra_true*ones(size(t)),'--k');
    plot(t,xhat(3,:));
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('Rv (mmHg/mL/s)');         
    ax = [ax, gca];    
    
    linkaxes(ax,'x');