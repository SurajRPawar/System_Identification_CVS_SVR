function signals = comp_data_generator(parameters, num_cycles, vad, heart_condition)
%% Function file to generate computational model data
%{
--------------------------- Description -----------------------------------
Simulate the 12 state computational model of the CVS (Gohean 2013) and
output generated signals. These signals to be used for system identification
tests. 

Notes
- A clean model is simulated without process noise. 
- Separately, another function will add additive white Gaussian
measurement noise. 
- For UKF, we will assume some values of process noise that will be tuned
heuristically.

----------------------------- Inputs --------------------------------------
1. Parameters    : (1) A       (mmHg)
                   (2) B       (1/mL)
                   (3) Emax    (mmHg/mL)
                   (4) Cs      (mL/mmHg)
                   (5) Rsvr    (mmHg/mL/s)
                   (6) Rv      (mmHg/mL/s)
                   (7) HR      (bpm)

2. num_cycles    : number of heart cycles to simulate
3. vad           : 1 for vad on, 0 for vad off                
5. heart_condition : 1 = healthy, 2 = heart failure. This is used to set
                     initial conditions
                 
----------------------------- Outputs -------------------------------------
signals : (1) Plv, Left ventricular pressure (mmHg)
          (2) Pao, Aortic pressure (mmHg)
          (3) Pla, Left Atrial pressure (mmHg)
          (4) Pra, Right Atrial pressure (mmHg)
          (5) Vlv, Left ventricle volume (mL)
          (6) Qa, Aortic flow rate (mL/s)
          (7) Qvad, VAD flow rate (mL/s)
          (8) Time vector (s)

---------------------------- Versions -------------------------------------
%{
v1 : Suraj R Pawar, 7-14-2020
    - Initialize
%}
v2 : Suraj R Pawar, 7-20-2020
    - Remove process noise based monte carlo simulation
    - Only perform clean Euler Simulation
    - Additive white Gaussian measurement noise will be added later in a
    separate function
%}

%% Input Parameters

    % Timing
    HR = parameters(7);
    tend = 60/HR*num_cycles;        
    dt = 0.0001;             % Simulation time step
    t = [0:dt:tend];
    rho = 1.055*1000;       %Density of blood (?) [kg/m^3]
    mu = 0.0027;          

    % Noise parameters (not used)
    %{
    pressure_noise_std = process_noise(1); % Pressure noise standard deviation
    volume_noise_std = process_noise(2);   % Volume noise standard deviation
    flow_noise_std = process_noise(3);     % Flow noise standard deviation    
    process_noise_std = [
                         volume_noise_std    % VLA
                         volume_noise_std    % VLV
                         pressure_noise_std  % PSA
                         flow_noise_std      % QSA
                         pressure_noise_std  % PST
                         pressure_noise_std; % PSV
                         volume_noise_std    % VRA
                         volume_noise_std    % VRV
                         pressure_noise_std  % PPA
                         flow_noise_std      % QPA
                         pressure_noise_std  % PPT
                         pressure_noise_std  % PPV
                        ];
                    
%     meas_noise_std = [
%                       5;  % PLV
%                       5;  % PAO
%                       1;  % Qa
%                      ];
                 
    %}
    
    % Initial Conditions
    %{
    CVS states: x = [VLA, VLV, PSA, QSA, PST, PSV, VRA, VRV, PPA, QPA, PPT, PPV];
    %}
    
    % Healthy
    x0_c1 = [
            29.9653
            94.4357
            107.5747
            53.8608
            99.7147
            8.5037
            53.0351
            115.2163
            11.1374
            81.5000
            5.5765
            2.0068
             ];
    
    % Heart Failure
    x0_c2 = [
            79.5771
            282.9557
            70.2232
            22.1339
            66.9627
            5.4734
            42.5886
            110.2630
            20.0657
            51.1726
            16.5074
            14.2768
            ];
    
    if heart_condition == 1 % Healthy
        x0 = x0_c1;  
    else                    % Heart Failure
        x0 = x0_c2;
    end

%% Steady state simulation

    t_short = [0:dt:60/HR*50]; % Reach steady state 
    
    %fprintf('Beginning simulation to find steady state \n'); tic;    
        
    % ODE simulation
    [T,X]   = ode15s(@(t, x)CVS_ODE(t,x,rho,parameters,vad),t_short,x0); 

    % Sort variables
    X_short = X;        
    x0_new = X_short(end,:)';            
    
    %fprintf('Finished steady state simulation in %.2f seconds\n', toc);

%% Euler Sim
   
    % Memory allocation
    steps = length(t);
    Xeu   = zeros(12, steps); 
    Xeu(:,1) = x0_new;    
    PLV = zeros(steps,1); 
    PLA = PLV;
    PSA = PLV;
    PRA = PLV;
    Qvad = zeros(steps, 1);
    QA = PLV;
    QM = PLV;    
    en = PLV;                    
    
    for i = 1:steps
        % State Dynamics propagation
        [xdot,y] = CVS_ODE(t(i),Xeu(:,i),rho,parameters,vad);
        if i ~= steps
            Xeu(:,i+1) = Xeu(:,i) + xdot.*dt;
        end

        % Collect outputs
        PSA(i) = y(1);
        PLV(i) = y(2);
        QA(i) = y(3);
        QM(i) = y(4);
        en(i) = y(5);
        Qvad(i) = y(6);
        PLA(i) = y(7);
        PRA(i) = y(8);
    end   
        
    VLV = Xeu(2,:).';    
    
    %{
    signals : (1) Plv, Left ventricular pressure (mmHg)
              (2) Pao, Aortic pressure (mmHg)
              (3) Pla, Left Atrial pressure (mmHg)
              (4) Pra, Right Atrial pressure (mmHg)
              (5) Vlv, Left ventricle volume (mL)
              (6) Qa, Aortic flow rate (mL/s)
              (7) Qvad, VAD flow rate (mL/s)
    %}
    signals = [PLV, PSA, PLA, PRA, VLV, QA, Qvad, t.'];
end