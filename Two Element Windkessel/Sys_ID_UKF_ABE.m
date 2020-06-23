%% A, B, Emax estimation using UKF and Two Element Windkessel model
%{
---------------------------- Description ----------------------------------
Uses Unscented Kalman Filter (UKF) to perform system identification of the
CVS using a two element model as the model basis. Sequential estimation
will be performed. First, we estimate B and Emax during ejection. Then, we
estimate A during filling. During each estimation phase, we use the last
best estimate of the other parameters. This repetitive process goes on for
all beats.

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
1. We assume that we know Pr, SVR, Cao, Rv
2. SVR can be estimated using the formula
        SVR = (MAP - CVP) / CO
3. We estimate the value of Rv from Pao differential equation

-------------------------------- Verisons ---------------------------------
%{
v1 : 6-10-2020, Suraj R Pawar
    - Initialize
v2 : Suraj R Pawar, 6-11-2020
    - Tested after adding 'counts' variable to func_handle_parameters
v3 : Suraj R Pawar, 6-12-2020
    - Added Qafilter structure containing filtered signal and thresholds
%}
v4 : Suraj R Pawar, 6-22-2020
    - Experiment with Healthy computational model data
    - Switched to a square root law for the aortic flow
%}

close all; clear all; clc;
include_us;

%% ----------------------- User Inputs ------------------------------------

    % Data file to load
    datafilename = 'sync_1_1_SysID.mat';
    
    % Set values for Cs, Rv and Pr
    Cs_true = 1.520; % (mL/mmHg)
    Pr_true = 14;    % (mmHg)
    Rv = 0.006;      % (mmHg/mL/s)
    
    % Set Qa filtering parameters and thresholds
    Qalowpassf = 8;    % Low pass cutoff frequency (Hz)
    Qaupper = 200;      % Where to set Qcheck flag to true
    Qalower = 0;        % When to stop marking ejection
    
    % Initial Conditions for B, E estimator
    V0_guess_be = 45;                       % Guess for unstressed blood volume
    Vlv0_be = 45;                         % Initial left ventricle volume (mL)
    Ps0_be = 65;                           % Initial systemic pressure (mmHg)
    B0 = 0.1;
    Emax0 = 10;
    
    % Initial conditions for A estimator
    V0_guess_a = 45;
    Vlv0_a = Vlv0_be;
    Ps0_a = 97;            
    A0 = 0.01;    
    
    % Initial covariances
    p0_states = diag([1, 5]);               % Initial error covariance of states Vbar and Ps
    p0_be = diag([1e-3, 5]);                % Initial error covariances for parameters B, Emax
    p0_a = [5];                             % Initial error covariance for parameter A
    
    % Process noise terms
    data.process_noise_std = [1; 5];        % Use only if you want to override default noise terms
    param_noise_std_be = [1e-10; 1e-10];    % White noise standard deviation for parameters [B, Emax]
    param_noise_std_a = [1e-10];            % White noise standard deviation for parameters [A]
    
    % UKF parameters (common for all parameter estimators)
    alpha = 1e-3;
    kappa = 0;
    beta = 2;                          % 2 is optimal for Gaussian distributions
    
    % What is the version of this experiment ? 
    experiment_versions;               % Contains description of versions  
    
    % Show GUI progress bar or not
    waitflag = 0;
    
%% ------------------ Measurements and parameters -------------------------
    data = noisy_twoelem_data(datafilename);      % Load noisy two element data using this function
      
    % True value of parameters
    
    Rsvr_true = data.parameters(2);    
    Ra_true = Rv;
    Rm_true = Rv;
    A_true = data.parameters(6);    
    B_true = data.parameters(7);
    E_true = data.parameters(8);
    V0_true = data.parameters(9);
    HR_true = data.parameters(10);
    t_c_true = data.parameters(11);
    t_vc_true = 0.6; %data.parameters(12);    
    parameters = [Rsvr_true; Cs_true; Pr_true; Ra_true; Rm_true; HR_true; t_vc_true; t_c_true;...
                  A_true; B_true; E_true; V0_true]; % To be passed to UKF       

    % Filter Qa
    dt_original = data.dt;
    Fs = 1/dt_original;                             % Sampling frequency (Hz)
    Qa_original = data.Qa;
    Qa_filtered = lowpass(Qa_original, Qalowpassf, Fs,...
                          'ImpulseResponse','iir','Steepness',0.9);
    
    % Interpolate all measurements to 1ms timing    
    dt = 0.001;
    tf = data.num_beats*data.parameters(11);        % Final time for simulation
    t_original = [0 : dt_original : tf];            % Time vector from measurement file
    t = [0: dt : tf];                               % Time vector with 1 ms time steps    
    Ps_original = data.Ps;    
    Plv_original = data.Plv;
    Qvad_original = data.Qvad;    
    Vbar_original = data.x(1,:) - V0_true;    
    
    Plv = interp1(t_original, Plv_original, t, 'linear', 'extrap');
    Pao = interp1(t_original, Ps_original, t, 'linear', 'extrap');
    Qa = interp1(t_original, Qa_original, t, 'linear', 'extrap');
    Qa_filtered = interp1(t_original, Qa_filtered, t, 'linear', 'extrap');
    Qvad = interp1(t_original, Qvad_original, t, 'linear', 'extrap');
    Vbar = interp1(t_original, Vbar_original, t, 'linear', 'extrap');   % Vbar = Vlv - V0
    
    y = [Plv; Pao; Qa];    
    
    Qafiltstruct.signal = Qa_filtered;
    Qafiltstruct.upper_threshold = Qaupper;
    Qafiltstruct.lower_threshold = Qalower;
    
    % Noise and parameters for UKF
    %{
    Further comments : 
    This is the white Gaussian noise strength for the continuous time
    model. We need to convert this into Ito's form and integrate using
    Euler integration. When we do convert it to Ito's form inside the
    function, the variance of the Brownian motion is q*dt
    %}
    
    q_be = diag([data.process_noise_std.^2; param_noise_std_be.^2]);    % Process noise covariance for B, E estimator
    q_a = diag([data.process_noise_std.^2; param_noise_std_a.^2]);      % Process noise covariance for A estimator
    
    r = diag(data.meas_noise_std.^2);                                   % Measurement noise covariance
    ukf_params = [alpha; kappa; beta];
    
    x0_be = [Vlv0_be - V0_guess_be; Ps0_be];                            % First state is Vbar = Vlv - V0
    x0_a = [Vlv0_a - V0_guess_a; Ps0_a];
    
    theta0_be = [B0; Emax0];
    theta0_a = [A0];
    
%% -------------------------------- UKF -----------------------------------
    fprintf('Beginning UKF estimation \n'); tic;    
    parametersl = parameters;
    tbe_ukf = [];
    ta_ukf = [];
    xbe_ukf = [];
    xa_ukf = [];
    ybe_ukf = [];
    ya_ukf = [];
    A_ukf = [];
    B_ukf = [];
    Emax_ukf = [];
    Qa_ukf_a = [];
    Qa_ukf_be = [];
    
    if waitflag == 1 
        f = waitbar(0,'UKF'); 
        console_freq = floor(data.num_beats/3);
    end
    
    for i = 1 : data.num_beats          
        
        % Console out
            if waitflag == 1
                if mod(i,console_freq) == 0                    
                    waitbar((i/data.num_beats),f,'UKF');
                end
            end
                
        % Collect all samples from current beat
            t_start = (i-1)*t_c_true;
            t_stop = i*t_c_true;
            
            index_start = find(t >= t_start, 1);
            index_stop = index_start + find(t(index_start:end) <= t_stop, 1 ,'last') -1;

            tl = t(index_start:index_stop);            
            yl = y(:,[index_start:index_stop]);
            Qvadl = Qvad(index_start:index_stop);
            Qa_filteredl = Qa_filtered(index_start:index_stop);
            Qafiltstruct.signal = Qa_filteredl;
            
        % B, E estimator                        
            if i ~= 1
                p_be = squeeze(Paug_be([3:4],[3:4],end));   % Last error covariance of parameters
                parametersl(9) = xhat_a(3,end);             % Last good estimate of A 
                theta0l = xhat_be([3:4],end);               % Last good estimate of B, Emax                
                p0l_states = squeeze(Paug_a([1:2],[1:2],end));                
                x0_be = xhat_a([1:2],end);
            else
                p_be = p0_be;
                parametersl(9) = A0;
                theta0l = theta0_be;
                p0l_states = p0_states;                
            end
            version = 1;
            p0l = blkdiag(p0l_states, p_be);                % Reset P for states
            [xhat_be, yhat_be, Paug_be, t_be] = func_TwoElem_SysID_UKF_BE_Ejection(tl, yl, Qvadl, x0_be, theta0l, p0l,...
                                                                       q_be, r, parametersl, ukf_params,...
                                                                       version, Qafiltstruct);
                                                               
        % A estimator
            if i~= 1
                p_a = squeeze(Paug_a(3,3,end));                
                theta0l = xhat_a(3,end);                   
            else
                p_a = p0_a;                
                theta0l = A0;
                
            end
            x0_a = xhat_be([1:2],end);
            p0l_states = squeeze(Paug_be([1:2],[1:2],end));            
            parametersl(10) = xhat_be(3,end);
            parametersl(11) = xhat_be(4,end);
            version = 2;
            p0l = blkdiag(p0l_states, p_a);
            [xhat_a, yhat_a, Paug_a, t_a] = func_TwoElem_SysID_UKF_A(tl, yl, Qvadl, x0_a, theta0l, p0l,...
                                                                   q_a, r, parametersl, ukf_params,...
                                                                   version, Qafiltstruct);
                                                               
        % Store in variables
            tbe_ukf = [tbe_ukf, t_be];
            ta_ukf = [ta_ukf, t_a];
            
            xbe_ukf = [xbe_ukf, xhat_be([1:2],:)];
            xa_ukf = [xa_ukf, xhat_a([1:2],:)];
            
            ybe_ukf = [ybe_ukf, yhat_be];
            ya_ukf = [ya_ukf, yhat_a];
                                  
            A_ukf = [A_ukf, xhat_a(3,:)];
            B_ukf = [B_ukf, xhat_be(3,:)];
            Emax_ukf = [Emax_ukf, xhat_be(4,:)];    
            
            Qa_ukf_a = [Qa_ukf_a, yhat_a(3,:)];
            Qa_ukf_be = [Qa_ukf_be, yhat_be(3,:)];
    end
    
    if waitflag == 1
        close(f);
    end
        
    fprintf('UKF estimation finished in %.2f seconds\n', toc);    
    
%% -------------------- Calculate accuracies ------------------------------
    B_final_est = mean(B_ukf(end-1000:end));
    A_final_est = mean(A_ukf(end-1000:end));
    Emax_final_est = mean(Emax_ukf(end-1000:end));
    
    B_acc = (1/B_true)*(B_true - B_final_est)*100;
    A_acc = (1/A_true)*(A_true - A_final_est)*100;
    Emax_acc = (1/E_true)*(E_true - Emax_final_est)*100;
    
    fprintf('==========================================\n');
    fprintf('ESTIMATES \n');
    fprintf('A: %.3f, accuracy: %.2f%% \n', A_final_est, A_acc);
    fprintf('B: %.3f, accuracy: %.2f%% \n', B_final_est, B_acc);
    fprintf('Emax: %.3f, accuracy: %.2f%% \n', Emax_final_est, Emax_acc);
    
%% --------------------------- Figures ------------------------------------
    figure;       
    default_marker = 3;
    ax = [];
    
    subplot(3,2,1); % Vbar
    hold on;
    plot(t,Vbar);
    plot(tbe_ukf,xbe_ukf(1,:),'r.', 'MarkerSize',default_marker);
    plot(ta_ukf,xa_ukf(1,:),'g.', 'MarkerSize',default_marker);   
    hold off;
    legend({'Measured','BE UKF', 'A UKF'},'Orientation','horizontal');
    title('Vbar (mL)');
    xlabel('Time (s)');
    ax = [ax, gca];
    
    subplot(3,2,2); % Plv
    hold on;
    plot(t,Plv);
    plot(tbe_ukf,ybe_ukf(1,:),'r.', 'MarkerSize',default_marker);
    plot(ta_ukf,ya_ukf(1,:),'g.', 'MarkerSize',default_marker);
    hold off;
    legend({'Measured','BE UKF', 'A UKF'},'Orientation','horizontal');
    title('Plv (mmHg)');
    xlabel('Time (s)');
    ax = [ax, gca];
    
    subplot(3,2,3); % Ps
    hold on;
    plot(t,Pao);
    plot(tbe_ukf,xbe_ukf(2,:),'r.', 'MarkerSize',default_marker);
    plot(ta_ukf,xa_ukf(2,:),'g.', 'MarkerSize',default_marker);
    hold off;
    legend({'Measured','BE UKF', 'A UKF'},'Orientation','horizontal');
    title('Pao (mmHg)');
    xlabel('Time (s)');    
    ax = [ax, gca];
    
    subplot(3,2,4); % A
    hold on;
    plot(t,A_true*ones(size(t)),'--k');
    plot(ta_ukf,A_ukf);
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('A (mmHg)');
    xlabel('Time (s)');       
    ax = [ax, gca];
    
    subplot(3,2,5); % B
    hold on;
    plot(t,B_true*ones(size(t)),'--k');
    plot(tbe_ukf,B_ukf);
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('B (1/mL)');
    xlabel('Time (s)');
    ax = [ax, gca];
    
    subplot(3,2,6); % Emax
    hold on;
    plot(t,E_true*ones(size(t)),'--k');
    plot(tbe_ukf,Emax_ukf);
    hold off;
    legend({'Actual','Estimated'},'Orientation','horizontal');
    title('Emax (mmHg/mL)');
    xlabel('Time (s)');
    ax = [ax, gca];
    
    linkaxes(ax, 'x');