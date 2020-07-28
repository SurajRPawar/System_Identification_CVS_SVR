function [A_final_est, B_final_est, Emax_final_est, y] = ABE_estimator(noisy_signals, parameters, Pr_true, Cs_true, Rv, Rsvr, t_vc_true, process_noise, meas_noise, Qa_filter, initial_guesses, num_cycles)
%% A, B, Emax estimator function file
%{
---------------------------- Description ----------------------------------
Uses Unscented Kalman Filter (UKF) to perform system identification of the
CVS using a two element model as the model basis. Sequential estimation
will be performed. First, we estimate B and Emax during ejection. Then, we
estimate A during filling. During each estimation phase, we use the last
best estimate of the other parameters. This repetitive process goes on for
all beats.

----------------------------- Inputs --------------------------------------
noisy_signals  : Signals from the computational model data generator, with
                additive Gaussian white noise added. 
                    (1) Plv, Left ventricular pressure (mmHg)
                    (2) Pao, Aortic pressure (mmHg)
                    (3) Pla, Left Atrial pressure (mmHg)
                    (4) Pra, Right Atrial pressure (mmHg)
                    (5) Vlv, Left ventricle volume (mL)
                    (6) Qa, Aortic flow rate (mL/s)
                    (7) Qvad, VAD flow rate (mL/s)
                    (8) Time vector (s)
parameters     : Parameters used for computational model data generation
                    (1) A       (mmHg)
                    (2) B       (1/mL)
                    (3) Emax    (mmHg/mL)
                    (4) Cs      (mL/mmHg)
                    (5) Rsvr    (mmHg/mL/s)
                    (6) Rv      (mmHg/mL/s)
                    (7) HR      (bpm)
Pr_true        : True value of Pr (approximated from Plv signal) (mmHg)
Cs_true        : Systemic compliance (mL/mmHg)
Rv             : Valve resistances (mmHg/mL/s)
Rsvr           : Systemic Vascular Resistance (mmHg/mL/s)
t_vc_true      : Ventricular contaction time (s)
process_noise  : Standard deviation for process noise of pressure, volume
                 and flow
meas_noise     : Standard deviation of measurement noise for pressure, flow,
                 volume
Qa_filter      : Structure with two elements : 
                 lowpass : Frequency for lowpass filter (Hz)
                 upper : Threshold above which ejection starts (mL)
                 lower : Threshold below which ejection stops (mL)
initial_guesses : Structure with following elements: 
                  Vbar0 (mL)
                  Ps0 (mmHg)
                  A0 (mmHg)
                  B0 (1/mL)
                  Emax0 (mmHg)
num_cycles      : Number of heart cycles in the dataset

-------------------------------- Outputs ----------------------------------
A_final_est    : Final estimated value of A (mmHg)
B_final_est    : Final estimated value of B (1/mL)
Emax_final_est : Final estimated value of Emax (mmHg/mL)
y : Structure containing
    (1) tbe : Time vector of B and E estimates
    (2) xbe : States (Vlv; Pao) from B and E estimation phases
    (3) b   : Time history of B estimates
    (4) e   : Time history of E estimates
    (5) ybe : Estimated outputs (Plv, Pao, Qa) from B and E estimation
              phase
    (6) ta : Time vector from A estimation
    (7) xa : States from A estimation phase
    (8) ya : Estimated outputs from A esitmation phase
    (9) a  : Time history of A estimates
    
-------------------------------- Verisons ---------------------------------
%{
v1 : Suraj R Pawar, 7-16-2020
    - Initialize
%}
v2 : Suraj R Pawar, 7-20-2020
    - Added y structure output
%}

%% ----------------------- Prepare parameters -----------------------------           
    
    % Set Qa filtering parameters and thresholds
    Qalowpassf = Qa_filter.lowpass; % Low pass cutoff frequency (Hz)
    Qaupper = Qa_filter.upper;      % Where to set Qcheck flag to true
    Qalower = Qa_filter.lower;      % When to stop marking ejection
    
    % Initial Conditions for B, E estimator    
    Vbar0_be = initial_guesses.Vbar0;
    Ps0_be = initial_guesses.Ps0;
    B0 = initial_guesses.B0;
    Emax0 = initial_guesses.Emax0;
    
    % Initial conditions for A estimator
    Vbar0_a = initial_guesses.Vbar0;
    Ps0_a = initial_guesses.Ps0;
    A0 = initial_guesses.A0;
    
    % Initial covariances
    p0_states = diag([1, 1]);
    p0_be = diag([1e-3, 1]);
    p0_a = [1];
    
    % Process noise terms    
    param_noise_std_be = [1e-10; 1e-10];    % White noise standard deviation for parameters [B, Emax]
    param_noise_std_a = [1e-10];            % White noise standard deviation for parameters [A]
    
    % UKF parameters (common for all parameter estimators)
    alpha = 1e-3;
    kappa = 0;
    beta = 2;       % 2 is optimal for Gaussian distributions        
    
    % Show GUI progress bar or not
    waitflag = 0;
    
%% ------------------ Measurements and parameters -------------------------          
    % True value of parameters
    
    Rsvr_true = Rsvr;    
    Ra_true = Rv;
    Rm_true = Rv;
    A_true = parameters(1);
    B_true = parameters(2);
    E_true = parameters(3);    
    HR_true = parameters(7);
    t_c_true = 60/HR_true;
    V0_true = 5;    
    %t_vc_true = 0.6; %data.parameters(12);    
    
    parameters = [Rsvr_true; Cs_true; Pr_true; Ra_true; Rm_true; HR_true; t_vc_true; t_c_true;...
                  A_true; B_true; E_true; V0_true]; % To be passed to UKF       

    % Filter Qa
    dt_original = noisy_signals(2,8) - noisy_signals(1,8); % Finding dt from time vector
    Fs = 1/dt_original;                                    % Sampling frequency (Hz)
    Qa_original = noisy_signals(:,6).';
    Qa_filtered = lowpass(Qa_original, Qalowpassf, Fs,...
                          'ImpulseResponse','iir','Steepness',0.9);
    
    % Interpolate all measurements to 1ms timing    
    dt = 0.001;
    tf = noisy_signals(end,8);                      % Final time for simulation
    t_original = noisy_signals(:,8).';              % Time vector from measurement file
    t = [0: dt : tf];                               % Time vector with 1 ms time steps    
    Ps_original = noisy_signals(:,2).';
    Plv_original = noisy_signals(:,1).';
    Qvad_original = noisy_signals(:,7).';
    Vbar_original = noisy_signals(:,5).' - V0_true;    
    
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
    
    q_be = diag([process_noise(2).^2; process_noise(1).^2; param_noise_std_be.^2]); % Process noise covariance for B, E estimator
    q_a = diag([process_noise(2).^2; process_noise(1).^2; param_noise_std_a.^2]);   % Process noise covariance for A estimator
    
    r = diag([meas_noise(1).^2; meas_noise(1).^2; meas_noise(2).^2]);               % Measurement noise covariance
    ukf_params = [alpha; kappa; beta];
    
    x0_be = [Vbar0_be; Ps0_be]; % First state is Vbar = Vlv - V0
    x0_a = [Vbar0_a; Ps0_a];
    
    theta0_be = [B0; Emax0];
    theta0_a = [A0];
    
%% -------------------------------- UKF -----------------------------------
    %fprintf('Beginning UKF estimation \n'); tic;    
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
        
    for i = 1 : num_cycles         
        
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
        
    %fprintf('UKF estimation finished in %.2f seconds\n', toc);    
    
%% ------------------------- Outputs --------------------------------------
    y = struct();
    y.tbe = tbe_ukf;
    y.xbe = xbe_ukf;    
    y.b = B_ukf;
    y.e = Emax_ukf;
    y.ybe = ybe_ukf;
    y.ta = ta_ukf;        
    y.xa = xa_ukf;
    y.ya = ya_ukf;  
    y.a = A_ukf;
    
%% -------------------- Calculate accuracies ------------------------------
    B_final_est = mean(B_ukf(end-1000:end));
    A_final_est = mean(A_ukf(end-1000:end));
    Emax_final_est = mean(Emax_ukf(end-1000:end));
    
    B_acc = (1/B_true)*(B_true - B_final_est)*100;
    A_acc = (1/A_true)*(A_true - A_final_est)*100;
    Emax_acc = (1/E_true)*(E_true - Emax_final_est)*100;
    
    %{
    fprintf('==========================================\n');
    fprintf('ESTIMATES \n');
    fprintf('A: %.3f, accuracy: %.2f%% \n', A_final_est, A_acc);
    fprintf('B: %.3f, accuracy: %.2f%% \n', B_final_est, B_acc);
    fprintf('Emax: %.3f, accuracy: %.2f%% \n', Emax_final_est, Emax_acc);
    %}
%% --------------------------- Figures ------------------------------------
%{
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
    %}
end
