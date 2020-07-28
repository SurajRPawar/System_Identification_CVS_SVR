% Yao et al 2003 based parameter identifiability
%{
------------------------ Description --------------------------------------
Use the two element windkessel model and perform parameter identifiability.
Based on Yao et al. 2003. PhD Research -> SVR Phase 2 -> 8

To run this file, you need the following :
1. Computation_Meas.mat
2. two_elem_
-------------------------- Versions ---------------------------------------
%{
v1 : Suraj R Pawar, 7-28-2020
    - Initialize
%}
%}

clear all; close all; clc;
include_us;

%% User Inputs        
    % Loading of measurement file
        file_to_use = 1;                        % 1 = Computatinoal meas, 2 = Baseline new
        if file_to_use == 1
            data = load('Computation_Meas.mat');
        elseif file_to_use == 2
            data = load('Baseline_new.mat');
        end
    
    % Parameter Values
        % Systemic Circulation
            Cao = 1.8;               % Systemic compliance (mL/mmHg)
            Rsvr = 0.975;            % Systemic Vascular Resistance (mmHg.s/mL)

        % Pulmonary Circulation
            Pr  = 3.7;               % Pulmonary circulation approximated as constant pressure source (mmHg)

        % Valves (Aortic and Mitral)
            Rv = 0.0025;             % Valve resistance, kept same for both valves (mmHg.s/mL)

        % Left Ventricle parameters
            A  = 0.03;               % Minimum left ventricle elastance (mmHg/mL)
            B = 0.05;
            E = 3.25;                % Left ventricle maximum elastance (mmHg/mL)

        % Initial Conditions
            Vbar0 = 96.4641;  
            Pao0 = 114.0695;

        % Simulation parameters
            dt = 0.001;             % Euler integration time step (s)

        % For sensitivity calculation
            decimate = 1;           % After how many simulation iterations must the sensitivity coefficient matrix be calculated
            num_params = 7;         % Number of parameters included in estimability analysis
            num_meas = 3;           % Number of measurement signals used               
            
% ------------------------- END OF USER INPUTS ----------------------------

%% INITIAL CONDITIONS, INTERPOLATION and TIME PARAMETERS ------------------
        x0 = [Vbar0; Pao0];
        
        if file_to_use == 1
            HR = 80;                                % Heart rate (bpm). Use this when using the computational measurements file
            t_vc = (550 - 1.75*HR)/1000;            % Ventricular contraction time (s) for computational measurements file
            t_c = 60/HR;                            % Total cardiac cycle time (s) for computational measurements file
            num_cycles = floor(data.teu(end)/t_c);  % Total number of heart cycles contained in the data
            HR = HR*ones(1,num_cycles);
            t_vc = t_vc*ones(1,num_cycles);         % Form a t_vc array
            t_c = t_c*ones(1, num_cycles);          % Form a t_c array
        elseif file_t_use == 2
            HR = data.HR;                           % Heart rate (bpm)
            HR_avg = data.HR_avg;                   % Average heart rate (bpm)
            num_cycles = length(HR);                % Number of heart cycles to simulate        
            t_vc = 0.6*ones(1, num_cycles);         % Ventricular contraction time (s)
            t_c = (60/HR_avg)*ones(1, num_cycles);  % Total cardiac cycle time (s)
        end
        
    % Qvad input signal from the measurement file
        if file_to_use == 1
            Z_Qvad = data.Z_Qvad;
            Z_QAO = data.Z_CO - data.Z_Qvad;
        elseif file_to_use == 2
            Qvad = data.Qvad_data;                    
        end
    
    % Collect parameters in a variable
        parameters  = [Cao, Rsvr, Pr, Rv, A, B, E, t_c(1), t_vc(1)].';           
    
%% SIMULATION -------------------------------------------------------------
    
    % Simulation Parameters
        if file_to_use == 1
            t = [data.teu(1) : dt : data.teu(end)];
        elseif file_to_use == 2
            t = data.t_data;
            dt = t(2) - t(1);
        end
        
        parameters = [parameters; dt];
        
    % Interpolate measurements if required
        if file_to_use == 1
            Qvad = interp1(data.teu,Z_Qvad,t);
            Qao = interp1(data.teu,Z_QAO,t);
            Plv = interp1(data.teu,data.Z_PLV,t);
            Pao = interp1(data.teu,data.Z_PSA,t);
        elseif file_to_use == 2
            Plv = data.Plv_data;
            Pao = data.Pao_data;
            Qao = data.Qao_data;
            Qvad = data.Qvad_data;
        end
        
    % Create simulation variables
        steps = length(t);
        Qasim = zeros(1,steps);                     % Simulated aortic flow
        Qmsim = zeros(1,steps);                     % Simulated mitral valve flow
        Plvsim = zeros(1,steps);                    % Simulated left ventricle pressure
        e_n_sim = zeros(1,steps);                   % Simulated values of normalized elastance
        dPsim = zeros(1,steps);                     % Simulated delta Pressures
        x = zeros(2, steps); x(:,1) = x0;           % State vector array
        stage = zeros(1,steps);                     % Marker to signify stage of cardiac cycle
        z_iters = floor(steps/decimate);            % Number of steps after decimation
        Z = zeros(num_meas*z_iters, num_params);    % Sensitivity coefficient matrix
        stage_decimated = zeros(1, z_iters);        % Stages after decimation
        x_decimated = zeros(2,z_iters);             % Decimated state vector
        
    % Simulation of model
        j = 1;                  % Current heart cycle
        decimate_counter = decimate;   % Counter for Z calculation while decemating the iterations
        decimate_index = 0;
        t_bc = t(1);            % Start time of this heart cycle (s)
        
        fprintf('Beginning Euler simulation \n'); tic;
        for i = 1:steps
            % Update the t_c, t_vc components of parameter vector
                parameters(8:9) = [t_c(j), t_vc(j)];
                
            % Get xdot and y at current time step
                [xdot, y, Zout] = two_elem_cvs_sensitivity(t(i), x(:,i), parameters, Qvad(i));                
                Qasim(i) = y(1);
                Qmsim(i) = y(2);
                Plvsim(i) = y(3);
                e_n_sim(i) = y(4);
                dPsim(i) = y(5);
                stage(i) = y(6);
               
            % Euler integration
                if i ~= steps
                    x(:,i+1) = x(:,i) + xdot*dt;
                end
                
            % Advance heart cycle if required
                if t(i) >= t_bc + t_c(j);   % If current time is >= Start time of this cycle + cardiac cycle time for this cycle
                    j = j + 1;                  % Then, increment heart cycle number index
                    t_bc = t(i);                % Set start time of heart cycle to current time value
                end
                
            % Sensitivity coefficient
                if decimate_counter >= decimate
                    decimate_index = decimate_index + 1;
                    stage_decimated(decimate_index) = stage(i);         % Store the current stage 
                    x_decimated(:,decimate_index) = x(:,i);
                    start_z_index = (decimate_index-1)*num_meas + 1;
                    stop_z_index = (decimate_index)*num_meas;
                    Z([start_z_index:stop_z_index],:) = Zout;
                    decimate_counter = 1;
                else
                    decimate_counter = decimate_counter + 1;
                end
        end
        fprintf('Finished Euler simulation in %.3f seconds \n', toc);

    % Separate sensitivity coefficient matrix into stages
        fprintf('Beginning separation of sensitivity coefficient matrix into stages \n', tic);
        
        % Get indexes for each stage in the decimated variables
        
            indexes_isovolumic = find(stage_decimated == 1);
            indexes_ejection = find(stage_decimated == 2);
            indexes_filling = find(stage_decimated == 3);
        
        % Form sensitivity coefficient matrix for each stage
            Z_ejection = Z(indexes_ejection,:);     % Sensitivity coefficient matrix during ejection
            Z_filling = Z(indexes_filling,:);       % Sensitivity coefficient matrix during filling
            Z_isovolumic = Z(indexes_isovolumic,:); % Sensitivity coefficient matrix during isovolumic stage
            
        fprintf('Finished preparation of sensitivity coefficient matrices for each stage in %.3f seconds \n', toc);
                    
%% ESTIMABILITY ANALYSIS --------------------------------------------------
    
    cutoff = 1   % When to stop the iterations
    
    % Ejection
        fprintf('\n Estimability Analysis for ejection stage \n');
        i = 1;
        while 1
            if i == 1
                Z_ejection_norms = vecnorm(Z_ejection);                         % Step 1
                [max_z, index(i)] = max(Z_ejection_norms);                      % Step 2
                fprintf('Estimable parameter index %d\n', index(i));    
                Xl_ejection = Z_ejection(:,index(i));                           % Step 3
            else
                Z_ejection_estimate = Xl_ejection*pinv(Xl_ejection)*Z_ejection; % Step 4
                Residual_ejection = Z_ejection - Z_ejection_estimate;           % Step 5
                Residual_ejection_norms = vecnorm(Residual_ejection);           % Step 6
                [max_residual, index(i)] = max(Residual_ejection_norms);        
                Xl_ejection = [Xl_ejection, Z_ejection(:,index(i))];            % Step 7
                fprintf('Estimable parameter index %d with residual %.3f\n', index(i), max_residual);
                if max_residual <= cutoff || i >= 7
                    fprintf('------------------------------------\n');
                    break;
                end
            end  

            i = i + 1;
        end
    
    % Filling
        fprintf('\n Estimability Analysis for Filling stage \n');
    
        %{
        Modify the Z matrix to remove the rows related to Qao. This
        measurement is 0 during the filing stage, and there is no point
        in keeping this row as a part of the sensitivity matrix
        %}
            Z_filling(3:3:end,:) = [];      % Delete every 3rd row, corresponding to Qao measurement
        
        i = 1;
        while 1
            if i == 1
                Z_filling_norms = vecnorm(Z_filling);                           % Step 1
                [max_z, index(i)] = max(Z_filling_norms);                       % Step 2
                fprintf('Estimable parameter index %d\n', index(i));    
                Xl_filling = Z_filling(:,index(i));                             % Step 3
            else
                Z_filling_estimate = Xl_filling*pinv(Xl_filling)*Z_filling;     % Step 4
                Residual_filling = Z_filling - Z_filling_estimate;              % Step 5
                Residual_filling_norms = vecnorm(Residual_filling);             % Step 6
                [max_residual, index(i)] = max(Residual_filling_norms);        
                Xl_filling = [Xl_filling, Z_filling(:,index(i))];               % Step 7
                fprintf('Estimable parameter index %d with residual %.3f\n', index(i), max_residual);
                if max_residual <= cutoff || i >= 7
                    fprintf('------------------------------------\n');
                    break;
                end
            end  

            i = i + 1;
        end
    
    % Isovolumic
        fprintf('\n Estimability Analysis for isovolumic stage \n');
            
        %{
        Modify the Z matrix to remove the rows related to Qao. This
        measurement is 0 during the isovolumic stage, and there is no point
        in keeping this row as a part of the sensitivity matrix
        %}
        
            Z_isovolumic(3:3:end,:) = [];      % Delete every 3rd row, corresponding to Qao measurement
        
        i = 1;
        while 1
            if i == 1
                Z_isovolumic_norms = vecnorm(Z_isovolumic);                             % Step 1
                [max_z, index(i)] = max(Z_isovolumic_norms);                            % Step 2
                fprintf('Estimable parameter index %d\n', index(i));    
                Xl_isovolumic = Z_isovolumic(:,index(i));                               % Step 3
            else
                Z_isovolumic_estimate = Xl_isovolumic*pinv(Xl_isovolumic)*Z_isovolumic; % Step 4
                Residual_isovolumic = Z_isovolumic - Z_isovolumic_estimate;             % Step 5
                Residual_isovolumic_norms = vecnorm(Residual_isovolumic);               % Step 6
                [max_residual, index(i)] = max(Residual_isovolumic_norms);        
                Xl_isovolumic = [Xl_isovolumic, Z_isovolumic(:,index(i))];              % Step 7
                fprintf('Estimable parameter index %d with residual %.3f\n', index(i), max_residual);
                if max_residual <= cutoff || i >= 7
                    fprintf('------------------------------------\n');
                    break;
                end
            end  

            i = i + 1;
        end
        
%% FIGURES ----------------------------------------------------------------
    figure(1); % Comparison of pressures and aortic flow
    subplot(2,1,1);
    hold on;
    plot(t,Plv,'r-',t,Plvsim,'r--');
    plot(t,Pao,'b-',t,x(2,:),'b--');
    hold off;
    ax1 = gca;
    ylabel('Pressure (mmHg)'); 
    legend('Plv Meas','Plv Sim', 'Pao Meas', 'Pao Sim');
    
    subplot(2,1,2);
    plot(t,Qao,'k-',t,Qasim,'--');
    ax2 = gca;
    ylabel('Flow (mL/s)'); xlabel('Time (s)');
    legend('Qao Meas','Qao Sim');
    
    linkaxes([ax1; ax2],'x');     