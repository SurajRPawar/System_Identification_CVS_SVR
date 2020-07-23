% Compare computational model data to two element model

close all; clear all; clc;
include_us;


%% Load signals generated from computational model
load('signals_healthy_vad.mat');

%% Parameters for Two Element Model
    % Parameters
    A = 0.0219;
    B = 0.0539;
    Emax = 1.9;
    Cs = 1.46;
    Rsvr = 0.9877;
    V0 = 5;
    Pr = 3;
    Rv = 0.005305;
    
    % Initial Conditions
    Vlv0 = 91.58 + V0;
    Ps0 = 105.6;
    
    % Timing
    HR = 80;
    num_cycles = 30;
    
%% Simulate Two Element Model

    tc = 60/HR;
    tf = num_cycles * tc;
    dt = 0.001;
    t = [0 : dt : tf];
    tvc = (550 - 1.75*HR)/1000;
    
    x0 = [Vlv0; Ps0];
    parameters = [Cs; Rsvr; Pr; Rv; Rv; A; B; Emax; V0; HR; tc; tvc];
    
    [T,X] = ode15s(@(t,x) two_elem_cvs(t,x,parameters),t,x0);
    
%% Get outputs
    
    % Memory Allocation
    steps = numel(t);
    Plv = zeros(1,steps);
    Qa = zeros(1,steps);
    
    for i = 1 : steps
        [~,y] = two_elem_cvs(t(i),X(i,:).',parameters);
        Qa(i) = y(1);
        Plv(i) = y(3);
    end
    
%% Compare
    ct = signals(:,8).';
    cPlv = signals(:,1).';
    cPao = signals(:,2).';
    cQa = signals(:,6).';
    
    Pao = X(:,2).';
    
%% Figures
    set(0, 'DefaultLineLineWidth',1);
    set(0, 'DefaultLegendFontWeight','bold');
    linewidth = 0.8;
    
    figure;
    Plvplot = subplot(2,1,1);
    Paoplot = subplot(2,1,2);
    %Qaplot = subplot(3,1,3);
    
    axes(Plvplot);
    hold on;
    plot(ct,cPlv,'-k');
    plot(t,Plv,':k');
    hold off;
    ax1 = gca;
    apply_axis_properties(ax1, linewidth);
    legend('Computational','Two Element');
    ylabel('\boldmath $P_{lv}$ \bf (mmHg)');
    
    axes(Paoplot);
    hold on;
    plot(ct,cPao,'-k');
    plot(t,Pao,':k');
    hold off;
    ax2 = gca;
    apply_axis_properties(ax2, linewidth);
    legend('Computational','Two Element');
    ylabel('\boldmath $P_{ao}$ \bf (mmHg)');
    xlabel('\bf Time (s)');
    
    linkaxes([ax1, ax2], 'x');
    
    figure;
    hold on;
    plot(ct,cQa,'-k');
    plot(t,Qa,':k');
    hold off;
    ax = gca;
    apply_axis_properties(ax, linewidth);
    legend('Computational','Two Element');
    ylabel('\boldmath $Q_{a}$ \bf mL/s');