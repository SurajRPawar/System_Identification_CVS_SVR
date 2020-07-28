function RMSE = Variation_Compare(parameters_est, parameters, signals, noisy_signals, Vbar0, Pr_true, num_cycles)
% Calculates RMSE between CVS model and measurements / simulations
%{
--------------------------- Description -----------------------------------
Accepts estimated parameters from system identification procedure, and
simulates the simplified two element CVS model with the identified
parameters. Compares the output of the identified model with the measured /
simulated outputs. Then calculates the RMSE errors between the outputs.

- Normalized RMSE used = RMSE *100/ range
- Normalized RMSE is used to gauge how 'close' is the simplified model to
measurements or simulations from computational model

--------------------------- Inputs ----------------------------------------
parameters_est : Parameters estimated using UKF based system identification
                (1) A       (mmHg)
                (2) B       (1/mL)
                (3) Emax    (mmHg/mL)
                (4) Cs      (mL/mmHg)
                (5) Rsvr    (mmHg/mL/s)
                (6) Rv      (mmHg/mL/s)
parameters     : Actual parameters. This vector is only used to extract HR
                (1) A       (mmHg)
                (2) B       (1/mL)
                (3) Emax    (mmHg/mL)
                (4) Cs      (mL/mmHg)
                (5) Rsvr    (mmHg/mL/s)
                (6) Rv      (mmHg/mL/s)
                (7) HR      (bpm)  <------- * We only need this *
signals        : Clean signals generated using computational model
                 simulation
                (1) Plv, Left ventricular pressure (mmHg)
                (2) Pao, Aortic pressure (mmHg)
                (3) Pla, Left Atrial pressure (mmHg)
                (4) Pra, Right Atrial pressure (mmHg)
                (5) Vlv, Left ventricle volume (mL)
                (6) Qa, Aortic flow rate (mL/s)
                (7) Qvad, VAD flow rate (mL/s)
                (8) Time vector (s)
noisy_signals  : Clean signals with additive measurement noise. Only need
                 Pao signal from this to generate initial Pao value
                (1) Plv, Left ventricular pressure (mmHg)
                (2) Pao, Aortic pressure (mmHg) <------ * Only need this
                (3) Pla, Left Atrial pressure (mmHg)
                (4) Pra, Right Atrial pressure (mmHg)
                (5) Vlv, Left ventricle volume (mL)
                (6) Qa, Aortic flow rate (mL/s)
                (7) Qvad, VAD flow rate (mL/s)
                (8) Time vector (s)
Vbar0          : Initial condition for Vbar0 state (mL)
Pr_true        : Value of constant pulmonary circulation pressure (mmHg)
num_cycles     : Number of heart cycles to simulate

------------------------ Outputs ------------------------------------------
RMSE           : Normalized root mean square errors of : 
                (1) Plv (% of range)
                (2) Pao (% of range)
                (3) dP  (% of range)
                (4) Qa  (% of range)

--------------------------- Versions --------------------------------------
%{
v1 : Suraj R Pawar, 7-27-2020
    - Initialize
%}
%}
%% Parameters for Two Element Model
    % Parameters
    A = parameters_est(1);
    B = parameters_est(2);
    Emax = parameters_est(3);
    Cs = parameters_est(4);
    Rsvr = parameters_est(5);
    V0 = 0;
    Pr = Pr_true;
    Rv = parameters_est(6);
    
    % Initial Conditions
    Vlv0 = Vbar0;
    Ps0 = noisy_signals(1,2);
    
    % Timing
    HR = parameters(7);
    %num_cycles = 30;
    
%% Simulate Two Element Model

    tc = 60/HR;
    tf = num_cycles * tc;
    dt = 0.0001;
    t = [0 : dt : tf];
    tvc = (550 - 1.75*HR)/1000;
    %tvc = 0.6;
    
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
    
    RMSE = [rmse_Plv; rmse_Pao; rmse_dP; rmse_Qa];
    %{
    fprintf('RMSE of PLV = %.3f (%% of mean) \n', rmse_Plv);
    fprintf('RMSE of Pao = %.3f (%% of mean) \n', rmse_Pao);
    fprintf('RMSE of dP = %.3f (%% of mean) \n', rmse_dP);
    fprintf('RMSE of Qa = %.3f (%% of mean) \n', rmse_Qa);
    %}
    
%% Figures
    %{
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
    ylabel('\boldmath $Q_{a}$ \bf mL/s');
    
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
    title(['\bf RMSE = \boldmath $',num2str(rmse_Plv,2),'$ \%']);
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
    title(['\bf RMSE = \boldmath $',num2str(rmse_Pao,2),'$ \%']);
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
    title(['\bf RMSE = \boldmath $',num2str(rmse_dP,2),'$ \%']);
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
    title(['\bf RMSE = \boldmath $',num2str(rmse_Qa,2),'$ \%']);
    ax4 = gca;
    %ylim([min(ylim) 250]);
    apply_axis_properties(ax4, linewidth);
    ylabel('\boldmath $Q_{a}$ \bf mL/s');
    xlabel('\bf Time (s)');
    %legend('Computational','Two Element');
    
    linkaxes([ax1,ax2,ax3,ax4],'x');
    %}