function [Csmean, Rvmean, y] = CsRv_estimator(noisy_signals, parameters, Pr_true, V0_true, t_vc_true, process_noise, meas_noise)
%% Cs and Rv estimator
%{
---------------------------- Description ----------------------------------
Uses Unscented Kalman Filter (UKF) to perform system identification of the
CVS using a two element model as the model basis. Assume that we know all
other parameters except for Cs and Rv. Estimation of these parameters can 
be performed during all stages of the cardiac cycle. 

Assumptions : 
1. We assume that we know Rsvr, A, B, Emax, Pr (can be approximated or 
   obtained from another estimator)
2. UKF can be used to estimation A, B, Emax
3. SVR can be estimated using the formula
        SVR = (MAP - CVP) / CO
4. Pr can be approximated from the Plv value during filling stage
5. We have assumed that the aortic and mitral valve resistances are the
same.

---------------------------- Inputs ---------------------------------------
noisy_signals : Signals from the computational model data generator, with
                additive Gaussian white noise added. 
                (1) Plv, Left ventricular pressure (mmHg)
                (2) Pao, Aortic pressure (mmHg)
                (3) Pla, Left Atrial pressure (mmHg)
                (4) Pra, Right Atrial pressure (mmHg)
                (5) Vlv, Left ventricle volume (mL)
                (6) Qa, Aortic flow rate (mL/s)
                (7) Qvad, VAD flow rate (mL/s)
                (8) Time vector (s)
parameters    : Parameters used for computational model data generation
                (1) A       (mmHg)
                (2) B       (1/mL)
                (3) Emax    (mmHg/mL)
                (4) Cs      (mL/mmHg)
                (5) Rsvr    (mmHg/mL/s)
                (6) Rv      (mmHg/mL/s)
                (7) HR      (bpm)
Pr_true       : True value of Pr (approximated from Plv signal) (mmHg)
V0_true       : Constant unstressed blood volume (mL)
                !! This parameter is not actually used anywhere !!
t_vc_true     : Ventricular contaction time (s)
process_noise : Standard deviation for process noise of pressure, volume
                and flow
meas_noise    : Standard deviation of measurement noise for pressure, flow,
                volume

------------------------------ Outputs ------------------------------------
Csmean : Mean esitmated value of Cs (mL/mmHg)
Rvmean : Mean estimated value of aortic and mitral valves (mmHg/mL/s)
y : Structure containing
    (1) t : Time vector (s)
    (2) xhat : Vector of augmented states [Pao; Cshat; Rvhat]

-------------------------------- Versions ---------------------------------
%{
v1 : 7-15-2020, Suraj R Pawar
    - Initialize
v2 : Suraj R Pawar, 7-16-2020
    - Added comments and descriptions
%}
v3 : Suraj R Pawar, 7-20-2020
    - Added structure output 'y'
%}

%% ---------------------- User Inputs -------------------------------------                
    
    % Initial Conditions    
    Ps0 = noisy_signals(1,2);           % Initial systemic pressure (mmHg)            
    Cs0 = 1;                            % Initial guess for systemic compliance (mL/mmHg)    
    Rv0 = 0.001;    
    
    % Initial covariances
    p0 = diag([1, 1, 1e-5]);           % Initial error covariance
        
    % Process noise terms
    %data.process_noise_std(2) = 1;      % If you want to override default value
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
    %data = noisy_twoelem_data(datafilename);      % Load noisy two element data using this function
      
    % True value of parameters
    Cs_true = parameters(4);
    Rsvr_true = parameters(5);
    Ra_true = parameters(6);
    Rm_true = parameters(6);
    A_true = parameters(1);
    B_true = parameters(2);
    E_true = parameters(3);
    %V0_true = data.parameters(9);
    HR_true = parameters(7);
    t_c_true = 60/HR_true;
    %t_vc_true = 0.6;%data.parameters(12);    
    parameters = [Rsvr_true; Cs_true; Pr_true; Ra_true; Rm_true; HR_true; t_vc_true; t_c_true;...
                  A_true; B_true; E_true; V0_true]; % To be passed to UKF       
    
    
    % Interpolate all measurements to 1us timing    
    dt = 0.001;
    %dt_original = noisy_signals(:,8).';
    tf = noisy_signals(end,8);   % Final time for simulation
    t_original = noisy_signals(:,8).';       % Time vector from measurement file
    t = [0: dt : tf];                          % Time vector with 1 ms time steps    
    Plv_original = noisy_signals(:,1).';
    Ps_original = noisy_signals(:,2).';
    Qvad_original = noisy_signals(:,7).';
    Qa_original = noisy_signals(:,6).';
         
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
    meas_noise_std = meas_noise(1:2);
    Ps_noise_std = process_noise(1);
    q = diag([Ps_noise_std.^2; param_noise_std.^2]);  % Process noise covariance
    r = diag(meas_noise_std.^2);                      % Measurement noise covariance
    ukf_params = [alpha; kappa; beta];
    x0 = [Ps0];                                       
    theta0 = [Cs0; Rv0];
    
%% ----------------------------- UKF --------------------------------------
    %fprintf('Beginning UKF estimation \n'); tic;    
    [xhat, yhat, Paug] = func_TwoElem_SysID_UKF_CsRv(t, y, u, x0, theta0, p0,...
                                                                   q, r, parameters, ukf_params,...
                                                                   version, waitflag);
    %fprintf('UKF estimation finished in %.2f seconds\n', toc);

%% -------------------------- Console output ------------------------------
    Csmean = mean(xhat(2,[end-50:end]));   
    Rvmean = mean(xhat(3,[end-50:end]));   
    %fprintf('Estimated Cs : %.3f mL/mmHg\n', Csmean);    
    %fprintf('Estimated Rv : %.3f mmHg/mL/s\n', Rvmean);   
    
%% ------------------------------- Outputs --------------------------------
    y = struct();
    y.t = t;
    y.x = xhat;
    
%% ------------------------------ Figures ---------------------------------
    %{
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
    %}
end