function noisy_signals = add_meas_noise(signals,meas_noise)
% Takes clean measurement signals and adds noise
%{
--------------------- Description -----------------------------------------
Used for regression tests with the UKF for system identification. Adds
additive white Gaussian noise to clean measurement signals.

- Noise is not added to time vector
- Although only PLV, PAO and QA signals are used, noise is added to all
other measurements. 

------------------------- Inputs ------------------------------------------
signals     : Clean signals of
                (1) Plv, Left ventricular pressure (mmHg)
                (2) Pao, Aortic pressure (mmHg)
                (3) Pla, Left Atrial pressure (mmHg)
                (4) Pra, Right Atrial pressure (mmHg)
                (5) Vlv, Left ventricle volume (mL)
                (6) Qa, Aortic flow rate (mL/s)
                (7) Qvad, VAD flow rate (mL/s)
                (8) Time vector (s)

meas_noise  : Standard deviation for measurement noise to be added. 
                (1) pressure meas noise std
                (2) flow meas noise std
                (1) volume meas noise std

---------------------- Versions -------------------------------------------
%{
v1 : Suraj R Pawar, 7-15-2020
    - Initialize
%}
%}

% Measurement noise standard deviations
    pressure_meas_noise = meas_noise(1);
    flow_meas_noise = meas_noise(2);
    volume_meas_noise = meas_noise(3);
    
% Add noise to signals
    noise_size = size(signals(:,1)); % Size of the noise vector
    Plv_noisy = signals(:,1) + pressure_meas_noise*randn(noise_size);
    Pao_noisy = signals(:,2) + pressure_meas_noise*randn(noise_size);
    Pla_noisy = signals(:,3) + pressure_meas_noise*randn(noise_size);
    Pra_noisy = signals(:,4) + pressure_meas_noise*randn(noise_size);
    Vlv_noisy = signals(:,5) + volume_meas_noise*randn(noise_size);
    Qa_noisy = signals(:,6) + flow_meas_noise*randn(noise_size);
    Qvad_noisy = signals(:,7) + flow_meas_noise*randn(noise_size);
    
% Collect into output
    noisy_signals = [Plv_noisy,...
                     Pao_noisy,...
                     Pla_noisy,...
                     Pra_noisy,...
                     Vlv_noisy,...
                     Qa_noisy,...
                     Qvad_noisy,...
                     signals(:,8)];
end

