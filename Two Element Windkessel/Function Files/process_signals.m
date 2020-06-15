%% Signal processing on Qa 
%{
------------------------- Description -------------------------------------
UKF code uses Qa (Aortic flow) signal for determining stage of cardiac
cycle. It helps to filter the signal before passing it on to the UKF. This
script reads the signals from a file, displays the power spectrum and
allows the user to test lowpass filters for a certain cutoff frequency. 

------------------------- Versions ----------------------------------------
v1 : Suraj R Pawar, 6-12-2020
    - Initialize
%}

clear all; close all; clc;

%% Load noisy measurements
data = noisy_twoelem_data('noisy_data_comp_hf.mat');
Plv = data.Plv;
Pao = data.Ps;
Qa = data.Qa;
Qvad = data.Qvad;
dt = data.dt;
HR = data.parameters(10);
num_beats = data.num_beats;
tc = data.parameters(11);
tf = tc*num_beats;
t = [0:dt:tf];

%% Power Spectrum
N = length(t);
if mod(N,2) ~= 0
    % Make samples of even length
    Plv = Plv(1:end-1);
    Pao = Pao(1:end-1);
    Qa = Qa(1:end-1);
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
wpass_plv = 12;
wpass_pao = 12;
wpass_qa = 30;

% Lowpass filtering
Plvf = lowpass(Plv, wpass_plv, Fs);
Paof = lowpass(Pao, wpass_pao, Fs);
Qaf = lowpass(Qa, wpass_qa, Fs);

% Smoothing
Plvs = smooth(Plvf,50);
Paos = smooth(Paof,50);
Qas = smooth(Qaf,5);

% Figure
figure;
subplot(4,1,1);
hold on;
plot(t,Plv,'color',0.6*ones(1,3));
plot(t,Plvs,'b-');
hold off;
title('Plv (mmHg)');

subplot(4,1,2);
hold on;
plot(t,Pao,'color',0.6*ones(1,3));
plot(t,Paos,'b-');
hold off;
title('Pao (mmHg)');

subplot(4,1,3);
hold on;
plot(t,Qvad(1:end-1),'color',0.6*ones(1,3));
plot(t,Qas,'b-');
hold off;
legend('Qvad', 'Qa smooth');
title('Qa (mL/s)');

subplot(4,1,4);
plot(t,Plvs,t,Paos);
title('Plv and Pao smoothed');

%% Determine average valve resistance
steps = length(Plvf);
j = 1;

for i = 1:steps
    dP = Plvs(i) - Paos(i);
    Ql = Qas(i);
    
    if dP >= 0 && Ql > 0;
        Rvest(j) = dP/Ql; 
        j = j + 1;
    end
end

Rvm = mean(Rvest)
figure;
hold on;
plot(Rvest,'k.');
plot(Rvm*ones(size(Rvest)),'r');
hold off;

