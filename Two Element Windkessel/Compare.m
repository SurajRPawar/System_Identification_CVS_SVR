% Compare computational model data to two element model

close all; clear all; clc;
include_us;


%% Load signals generated from computational model
load('Results/Resulting Signals/Resulting_Signals_Animal_Sim.mat');

%% Parameters for Two Element Model
    % Parameters
%     A = 0.207;
%     B = 0.016;
%     Emax = 0.3;
%     Cs = 1.3;
%     Rsvr = 1.0752;
%     V0 = 0;
%     Pr = 14;
%     Rv = 0.01476;

     % Healthy (new model)
%     A = 0.0252;
%     B = 0.0492;
%     Emax = 3.2;
%     Cs = 1.31;
%     Rsvr = 0.9748;
%     V0 = 5;
%     Pr = 3;
%     Rv = 0.004786;
    
%     % HF (new model)
%     A = 0.501;
%     B = 0.0145;
%     Emax = 0.35;
%     Cs = 0.889;
%     Rsvr = 1.0752;
%     V0 = 5;
%     Pr = 13.8;
%     Rv = 0.01475;
    
    % Animal (new model)
    A = 0.741;
    B = 0.0314;
    Emax = 2.6;
    Cs = 1.26;
    Rsvr = 0.7684;
    V0 = 5;
    Pr = 18;
    Rv = 0.005496;
    
    % Initial Conditions
%     Vlv0 = 259.5 + V0;
%     Ps0 = 73.67;
    
%     % Healthy (new model)
%     Vlv0 = 97.22 + V0;
%     Ps0 = 108.2;
% 
%     % HF (new model)
%     Vlv0 = 225.6 + V0;
%     Ps0 = 72.95;
    
    % Animal (new model)
    Vlv0 = 103 + V0;
    Ps0 = 62;
    
    % Timing
%     HR = 90;
%     num_cycles = 30;
    
    HR = 70.6057;
    num_cycles = 50;
    
%% Simulate Two Element Model

    tc = 60/HR;
    tf = num_cycles * tc;
    dt = 0.001;
    t = [0 : dt : tf];
    tvc = (550 - 1.75*HR)/1000;
    tvc = 0.6;
    
    x0 = [Vlv0; Ps0];
    parameters = [Cs; Rsvr; Pr; Rv; Rv; A; B; Emax; V0; HR; tc; tvc];
    Qvad_meas = signals(:,7).';
    t_meas = signals(:,8).';
    
    [T,X] = ode23s(@(t,x) two_elem_cvs(t,x,parameters, t_meas, Qvad_meas),t,x0);
    
%% Get outputs
    
    % Memory Allocation
    steps = numel(t);
    Plv = zeros(1,steps);
    Qa = zeros(1,steps);
    
    for i = 1 : steps
        [~,y] = two_elem_cvs(t(i),X(i,:).',parameters, t_meas, Qvad_meas);
        Qa(i) = y(1);
        Plv(i) = y(3);
    end
    
%% Compare
    ct = signals(:,8).';
    cPlv = signals(:,1).';
    cPao = signals(:,2).';
    cQa = signals(:,6).';
    
    Pao = X(:,2).';
    
%% RMSEs
    dPhat = Plv - Pao;
    dP = cPlv - cPao;
    se_dP = (dP - dPhat).^2;
    mse_dP = mean(se_dP);
    rmse_dP = sqrt(mse_dP);
    den_dP = range(dP);
    rmse_dP = rmse_dP*100/den_dP;
    
    se_Plv = (cPlv - Plv).^2;
    mse_Plv = mean(se_Plv);
    rmse_Plv = sqrt(mse_Plv);
    den_Plv = range(cPlv);
    rmse_Plv = rmse_Plv*100/den_Plv;
    
    se_Pao = (cPao - Pao).^2;
    mse_Pao = mean(se_Pao);
    rmse_Pao = sqrt(mse_Pao);
    den_Pao = range(cPao);
    rmse_Pao = rmse_Pao*100/den_Pao;
    
    se_Qa = (cQa - Qa).^2;
    mse_Qa = mean(se_Qa);
    rmse_Qa = sqrt(mse_Qa);
    den_Qa = range(cQa);
    rmse_Qa = rmse_Qa*100/den_Qa;
    
    fprintf('nRMSE of PLV = %.3f (%% of mean) \n', rmse_Plv);
    fprintf('nRMSE of Pao = %.3f (%% of mean) \n', rmse_Pao);
    fprintf('nRMSE of dP = %.3f (%% of mean) \n', rmse_dP);
    fprintf('nRMSE of Qa = %.3f (%% of mean) \n', rmse_Qa);
    
%% Figures
    set(0, 'DefaultLineLineWidth',1);
    set(0, 'DefaultLegendFontWeight','bold');
    set(0, 'DefaultLegendOrientation','vertical');
    set(0, 'DefaultLegendBox','off');
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
    plot(t,Qa,':r');
    hold off;
    ax = gca;
    apply_axis_properties(ax, linewidth);
    legend('Computational','Two Element');
    ylabel('\boldmath $Q_{a}$ \bf (mL/s)');
    
    figure; % Comparison of outputs
    decimation = 100;
    Plv_compare = subplot(2,2,1);
    Pao_compare = subplot(2,2,2);
    dP_compare = subplot(2,2,3);
    Qa_compare = subplot(2,2,4);
    
    axes(Plv_compare);
    hold on;
    plot(t,cPlv,'-k','LineWidth',2,'Color',[0.7,0.7,0.7]);
    plot(t(1:decimation:end),Plv(1:decimation:end),':k');
    hold off;
    title(['\bf nRMSE = \boldmath $',num2str(rmse_Plv,2),'$ \%']);
    ax1 = gca;
    %ylim([min(ylim) 100]);
    apply_axis_properties(ax1, linewidth);
    ylabel('\boldmath $P_{lv}$ \bf (mmHg)');
    legend('Computational','Two Element');
    
    axes(Pao_compare);
    hold on;
    plot(t,cPao,'-k','LineWidth',2,'Color',[0.7,0.7,0.7]);
    plot(t(1:decimation:end),Pao(1:decimation:end),':k');
    hold off;
    title(['\bf nRMSE = \boldmath $',num2str(rmse_Pao,2),'$ \%']);
    ax2 = gca;
    %ylim([min(ylim) 110]);
    apply_axis_properties(ax2, linewidth);
    ylabel('\boldmath $P_{ao}$ \bf (mmHg)');
    %legend('Computational','Two Element');
    
    axes(dP_compare);
    hold on;
    plot(t,dP,'-k','LineWidth',2,'Color',[0.7,0.7,0.7]);
    plot(t(1:decimation:end),dPhat(1:decimation:end),':k');
    hold off;
    title(['\bf nRMSE = \boldmath $',num2str(rmse_dP,2),'$ \%']);
    ax3 = gca;
    %ylim([min(ylim) 50]);
    apply_axis_properties(ax3, linewidth);
    ylabel('\boldmath $\Delta P$ \bf (mmHg)');
    xlabel('\bf Time (s)');
    %legend('Computational','Two Element');
    
    axes(Qa_compare);
    hold on;
    plot(t,cQa,'-k','LineWidth',2,'Color',[0.7,0.7,0.7]);
    plot(t(1:decimation:end),Qa(1:decimation:end),':k');
    hold off;
    title(['\bf nRMSE = \boldmath $',num2str(rmse_Qa,2),'$ \%']);
    ax4 = gca;
    %ylim([min(ylim) 250]);
    apply_axis_properties(ax4, linewidth);
    ylabel('\boldmath $Q_{a}$ \bf (mL/s)');
    xlabel('\bf Time (s)');
    %legend('Computational','Two Element');
    
    linkaxes([ax1,ax2,ax3,ax4],'x');