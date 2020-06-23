%% Signal processing on Qa 
%{
------------------------- Description -------------------------------------
UKF code uses Qa (Aortic flow) signal for determining stage of cardiac
cycle. It helps to filter the signal before passing it on to the UKF. This
script reads the signals from a file, displays the power spectrum and
allows the user to test lowpass filters for a certain cutoff frequency. 

------------------------- Versions ----------------------------------------
%{
v1 : Suraj R Pawar, 6-12-2020
    - Initialize
v2 : Suraj R Pawar, 6-16-2020
    - Added include_us script
%}
v3 : Suraj R Pawar, 6-23-2020
    - Signal trimming is added within the noisy_twoelem_data funciton file,
    therefore we don't need to do any signal trimming in this script. 
    - Removed the code section related to valve resistance calculation.
    That value is now estimated using the Cs equations
%}

clear all; close all; clc;

%% Load noisy measurements
data = noisy_twoelem_data('sync_1_1_SysID.mat');
Plv = data.Plv;
Pao = data.Ps;
Qa = data.Qa;
Qvad = data.Qvad;
dt = data.dt;
HR = data.parameters(10);
num_beats = data.num_beats;
tc = data.parameters(11);
tf = tc*num_beats;
t = data.t;

%% Power Spectrum
N = length(t);
if mod(N,2) ~= 0
    % Make samples of even length
    Plv = Plv(1:end-1);
    Pao = Pao(1:end-1);
    Qa = Qa(1:end-1);
    Qvad = Qvad(1:end-1);
    t = t(1:end-1);
    N = N-1;
end

Fs = 1/dt;

y_plv = fft(Plv); 
y_plv = y_plv(1:N/2 + 1);

y_pao = fft(Pao);
y_pao = y_pao(1:N/2 + 1);

y_qa = fft(Qa);
y_qa = y_qa(1:N/2 + 1);

% Power Spectrum
p_plv = abs(y_plv).^2/(Fs*N);
p_plv(2:end-1) = 2*p_plv(2:end-1);

p_pao = abs(y_pao).^2/(Fs*N);
p_pao(2:end-1) = 2*p_pao(2:end-1);

p_qa = abs(y_qa).^2/(Fs*N);
p_qa(2:end-1) = 2*p_qa(2:end-1);

% Figure
freq = [0 : Fs/N : Fs/2];
figure;
subplot(3,1,1);
plot(freq, 10*log10(p_plv));
title('PSD of PLV (dB/Hz)');
xlim([0 50]);

subplot(3,1,2);
plot(freq, 10*log10(p_pao));
title('PSD of PAO (dB/Hz)');
xlim([0 50]);

subplot(3,1,3);
plot(freq, 10*log10(p_qa));
title('PSD of QA (dB/Hz)');
xlabel('Frequency (Hz)');
xlim([0 50]);

%% Lowpass filtering

% Cutoff frequencies in Hz
wpass_plv = 8;
wpass_pao = 8;
wpass_qa = 8;

% Lowpass filtering
Plvf = lowpass(Plv, wpass_plv, Fs, 'ImpulseResponse','iir','Steepness',0.8);
Paof = lowpass(Pao, wpass_pao, Fs, 'ImpulseResponse','iir','Steepness',0.8);
Qaf = lowpass(Qa, wpass_qa, Fs, 'ImpulseResponse','iir','Steepness',0.98);

% Smoothing
Plvs = smooth(Plvf,10);
Paos = smooth(Paof,10);
Qas = smooth(Qaf,10);

% Figure
figure;
linewidth = 2;
ax = [];

subplot(4,1,1);
hold on;
plot(t,Plv,'color',0.7*ones(1,3), 'LineWidth',linewidth);
plot(t,Plvf,'b-');
hold off;
legend('Measured','Filtered');
title('Plv (mmHg)');
ax = [ax, gca];

subplot(4,1,2);
hold on;
plot(t,Pao,'color',0.7*ones(1,3), 'LineWidth',linewidth);
plot(t,Paof,'b-');
hold off;
legend('Measured','Filtered');
title('Pao (mmHg)');
ax = [ax, gca];

subplot(4,1,3);
hold on;
plot(t,Qvad,'k');
plot(t,Qa,'color',0.7*ones(1,3), 'LineWidth',linewidth);
plot(t,Qaf,'b-');
hold off;
legend('Qvad measured', 'Qa measured', 'Qa filtered');
title('Qa (mL/s)');
ax = [ax, gca];

subplot(4,1,4);
plot(t,Plvf,t,Paof);
title('Plv and Pao filtered');
ax = [ax, gca];

linkaxes(ax,'x');