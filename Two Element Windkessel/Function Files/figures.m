function [] = figures(signals, noisy_signals, y1, y2, V0_true, decimation, figurestyle, means)
% Figures
%{
------------------------------ Description --------------------------------
Generate figures for publication purposes. This function needs to be used
along with System Identification Regression tests using UKF based on Two
Element Windkessel model. 

---------------------------- Inputs ---------------------------------------
signals      : Clean signals generated from computational model data
                  (1) Plv, Left ventricular pressure (mmHg)
                  (2) Pao, Aortic pressure (mmHg)
                  (3) Pla, Left Atrial pressure (mmHg)
                  (4) Pra, Right Atrial pressure (mmHg)
                  (5) Vlv, Left ventricle volume (mL)
                  (6) Qa, Aortic flow rate (mL/s)
                  (7) Qvad, VAD flow rate (mL/s)
                  (8) Time vector (s)
noisy_signals : Same as signals, but with White Gaussian additive noise
y1 : Structure output from Cs and Rv estimator containing
    (1) t : Time vector (s)
    (2) xhat : Vector of augmented states [Pao; Cshat; Rvhat]
y2 : Structure output from A,B and E estimator containing
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
V0_true : Unstressed blood volume in LV (mL)
decimation : Decimation for plotting points
figurestyle : Style of figures
                1 = Individual figure windows for all signals
                2 - Predetermined grouping into subplots
means : Mean value of final estimates. This information is used to label
        the plots

--------------------------- Outputs ---------------------------------------
None

------------------------- Versions ----------------------------------------
%{
v1 : Suraj R Pawar, 7-20-2020
    - Initialize
%}
v3 : Suraj R Pawar, 7-20-2020
    - Added means as an input to label plots appropriately. 
%}

tc = y1.t;
xc = y1.x;

ta = y2.ta(1 : decimation : end);
xa = y2.xa(:, [1 : decimation : end]);
ya = y2.ya(:, [1 : decimation : end]);
a = y2.a(1 : decimation : end);

tbe = y2.tbe(1 : decimation : end);
xbe = y2.xbe(:, [1 : decimation : end]);
ybe = y2.ybe(:, [1 : decimation : end]);
b = y2.b(1 : decimation : end);
e = y2.e(1 : decimation : end);

figure;
subplot(2,2,1); % Vbar
hold on;
plot(noisy_signals(:,8), signals(:,5) - V0_true,'-k');
plot(ta, xa(1,:),'.r');
plot(tbe, xbe(1,:),'.r');
hold off;
title('Vbar (mL)');

subplot(2,2,2); % Plv
hold on;
plot(noisy_signals(:,8), signals(:,1),'-k');
plot(ta, ya(1,:),'.r');
plot(tbe, ybe(1,:),'.r');
hold off;
title('Plv (mmHg)');

subplot(2,2,3); % Pao
hold on;
plot(noisy_signals(:,8), signals(:,2),'-k');
plot(ta, xa(2,:),'.r');
plot(tbe, xbe(2,:),'.r');
hold off;
title('Pao (mmHg)');

subplot(2,2,4); % Qao
hold on;
plot(noisy_signals(:,8), signals(:,6),'-k');
plot(ta, ya(3,:),'.r');
plot(tbe, ybe(3,:),'.r');
hold off;
title('Qa (mL/s)');

%% Parameter Estimates

    % Common line specifications
    linewidth = 1;
    fontsize = 12;
    axwidth = 0.7;
    pos = [259.3333 , 115.3333 , 685.0000 , 514.3333];

    if figurestyle == 2
        figure('Position', pos);
        Csplot = subplot(3,2,1);
        Rvplot = subplot(3,2,2);
        Aplot = subplot(3,2,3);
        Bplot = subplot(3,2,4);
        Eplot = subplot(3,2,5);
    end
    
    % Cs    
    if figurestyle == 2 
        axes(Csplot); 
    else
        figure;
    end
    plot(tc,xc(2,:),'k','LineWidth', linewidth); grid on;
    if figurestyle == 2
       title(['\boldmath $\hat{C}_s \rightarrow', num2str(means(4),3),'$ \bf mL/mmHg']); 
    else
    xlabel('\bf Time (s)'); ylabel('\boldmath $\hat{C}_s$ \bf (mL/mmHg)');    
    end   
    ax = gca;
    box on;
    ax.LineWidth = axwidth;
    ax.FontWeight = 'bold';    

    % Rv
    if figurestyle == 2 
        axes(Rvplot);
    else
        figure; 
    end    
    plot(tc,xc(3,:),'k','LineWidth',linewidth); grid on;
    if figurestyle == 2
        title(['\boldmath $\hat{R}_v \rightarrow', num2str(means(6),4), '$ \bf mmHg.s/mL']);
    else
        xlabel('\bf Time (s)'); ylabel('\boldmath $\hat{R}_v$ \bf (mmHg.s/mL)');
    end    
    ax = gca;
    box on;
    ax.LineWidth = axwidth;
    ax.FontWeight = 'bold';    

  	% A
    if figurestyle == 2
        axes(Aplot);
    else
        figure;
    end
    plot(ta,a,'k','LineWidth',linewidth); grid on;
    if figurestyle == 2
        title(['\boldmath $\hat{A} \rightarrow', num2str(means(1),3) ,'$ \bf mmHg']);        
    else
        xlabel('\bf Time (s)'); ylabel('\boldmath $\hat{A}$ \bf(mmHg)');
    end
    ax = gca;
    box on;
    ax.LineWidth = axwidth;
    ax.FontWeight = 'bold';

    % B
    if figurestyle == 2
        axes(Bplot);
    else
        figure;
    end
    plot(tbe,b,'k','LineWidth',linewidth); grid on;
    if figurestyle == 2
        title(['\boldmath $\hat{B} \rightarrow',num2str(means(2),3),'$ \bf 1/mL']);
        xlabel('\bf Time (s)');
    else
        xlabel('\bf Time (s)'); ylabel('\boldmath $\hat{B}$ \bf (1/mL)');
    end
    ax = gca;
    box on;
    ax.LineWidth = axwidth;
    ax.FontWeight = 'bold';

    % Emax
    if figurestyle == 2
        axes(Eplot);
    else
        figure;
    end
    plot(tbe,e,'k','LineWidth',linewidth); grid on;
    if figurestyle == 2
        title(['\boldmath $ \hat{E}_{max} \rightarrow', num2str(means(3),2), '$ \bf mmHg/mL ']);
        xlabel('\bf Time (s)');
    else        
        xlabel('\bf Time (s)'); ylabel('\boldmath $\hat{E}_{max}$ \bf (mmHg)');
    end
    ax = gca;
    box on;
    ax.LineWidth = axwidth;
    ax.FontWeight = 'bold';
    end

