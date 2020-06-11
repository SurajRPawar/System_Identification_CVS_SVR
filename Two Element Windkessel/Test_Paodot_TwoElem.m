%% Test the Paodot state equation through all stages of CVS
%{
----------------------- Description ---------------------------------------
The Psdot state equation for the Two Element Windkessel model is the same
for all stages - filling, ejection and isovolumic expansion. Before I can
write a UKF code with that state equation for the estimation of Pr and Cs,
I wanted to make sure that we get the same Ps signal using this state
equation to perorm Euler Integration. 

Conclusion - Yes, this state equation is able to replicate the Ps
measured signal.

------------------------- Versions ----------------------------------------
v1 : Suraj R Pawar, 6-11-2020
    - Initialize
%}

clear all; close all; clc;
include_us;

%% Load Measurements
data = noisy_twoelem_data;
t0 = 0;
dt = data.dt;
tc = data.parameters(11);
num_beats = data.num_beats;
tf = tc * num_beats;
t = [t0 : dt : tf];

Rsvr = 1*data.parameters(2);
Cs = 1*data.parameters(1);
Pr = 1*data.parameters(3);

Plv = data.Plv;     % (mmHg)
Ps = data.Ps;       % (mmHg)
Qa = data.Qa;       % (mL/s)
Qvad = data.Qvad;   % (mL/s)
Qao = Qa + Qvad;

%% Euler Integration of Psdot
steps = length(t);
Pshat = zeros(1,steps);
Pshat(1) = Ps(1);

for i = 2 : steps
    Psp = Ps(i-1);
    tp = t(i-1);
    Qvadp = Qvad(i-1);
    Qap = Qa(i-1);
    
    Psdot = (Pr-Psp)/(Rsvr*Cs) + (Qvadp + Qap)/Cs;
    Psc = Psp + Psdot*dt;
    
    Pshat(i) = Psc;
end

%% Figures

figure;
pressures = subplot(3,1,1);
flows = subplot(3,1,2);
Psplot = subplot(3,1,3);
ax = [];

axes(pressures);
hold on;
plot(t, Plv, 'r-');
plot(t, Ps, 'b-');
hold off;
ylabel('Pressure (mmHg)'); xlabel('Time (s)');
legend('Plv','Ps');
ax = [ax, gca];

axes(flows);
hold on;
plot(t, Qa, 'r-');
plot(t, Qvad, 'b-');
hold off;
ylabel('Flow (mL/s)'); xlabel('Time (s)');
legend('Qa','Qvad');
ax = [ax, gca];

axes(Psplot);
hold on;
plot(t(1:10:end), Ps(1:10:end), 'k-', 'LineWidth',1);
plot(t(1:10:end), Pshat(1:10:end),'r--','LineWidth',1);
hold off;
ylabel('Pressure (mmHg)'); xlabel('Time (s)');
legend('Ps hat', 'Ps meas');
ax = [ax, gca];

linkaxes(ax, 'x');