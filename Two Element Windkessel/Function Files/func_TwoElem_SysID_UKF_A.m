function [xhat, yhat, p_aug, tselect] = func_TwoElem_SysID_UKF_A(t, y, Qvad, x0, theta0,...
                                                                              p0, q, r, parameters,...
                                                                              ukf_params, version, Qafiltstruct, waitflag)        
% UKF for A estimation using Two Element Windkessel model
%{
----------------------------- DESCRIPTION ---------------------------------
Take measurement signals, detect isovolumic expansion and filling and 
apply UKF only to these windows. At the beginning of each window, we reset
the states to their initial value, and reset the covariances corresponsing
to the states to the initial covariance. By states, we refer to [Vbar; Ps].
The estimates and covariances of the parameters being estimated will carry
forward.

----------------------------- ASSUMPTIONS ---------------------------------
- We know SVR from the equation (MAP - CVP)/CO
- We know Cs, Pr, B, Emax, Rv
- Only unknown parameters to estimate are A

----------------------------- INPUTS --------------------------------------
t          : Time vector in the format [t0 : dt: tf]
y          : Measurement signals of Plv, Pao, Qa (aortic flow without VAD included)
Qvad       : VAD flow rate input (mL/s)
x0         : Initial guess for states
theta0     : initial guess for parameters [A; B; Emax; Rv]
p0         : Initial error covariance for states and parameters
q          : Process noise covariance (states and parameters)
r          : Measurement noise covariance
parameters : True value of all parameters
             1. Rsvr : Systemic Vascular Resistance (mmHg/mL/s)
             2. Cs   : Systemic compliance (mL/mmHg)
             3. Pr   : Constant pulmonary circulation pressure (mmHg)
             4. Ra   : Aortic valve resistance (mmHg/mL/s)
             5. Rm   : Mitral valve resistance (mmHg/mL/s)
             6. HR   : Heart rate (bpm)
             7. t_vc : Ventricular contraction time (s)
             8. t_c  : Cardiac cycle time (s)
             9. A    : Parameter for PLV generation (mmHg)
             10. B   : Parameter for PLV generation (1/mL)
             11. E   : Maximum LV elastance (mmHg/mL)
             12. V0  : Unstressed blood volume (mL)
ukf_params : [alpha; kappa; beta]
version    : The version of experiment to run. This verison will change the
             way parameters and their dynamic equations are set. Refer to
             the MATLAB file 'experiment_versions.m'
Qafiltstruct : Structure containing : 
                1. signal : Filtered Qa signal (mL/s)
                2. upper_threshold : Upper threshold for setting Qcheck
                                     flag to TRUE
                3. lower_threshold : Lower threshold for marking end of
                                     filling stage.
waitflag (o) : Wait flag (0 or 1) to display GUI waitbar. 1 = display
               waitbar

----------------------------- OUTPUTS -------------------------------------
xhat    : Estimated states
theta   : Estimated parameters
p_aug   : Covariances for the augmented state [x; theta; v; n]
tselect : Time vector for iso contraction and ejection windows

----------------------------- VERSION -------------------------------------
%{
v1 : Suraj R Pawar, 5-26-2020
    - Initialize
v2 : Suraj R Pawar, 5-28-2020
    - Clean name of function
    - Clean documentation and comments
v3 : Suraj R Pawar, 5-28-2020
    - Add version mechanism to set parameter values and dynamic equations
    based on experiment being run
v4 : Suraj R Pawar, 6-9-2020
    - Fixed comments
v5 : Suraj R Pawar, 6-10-2020
    - Using Qa filtered signal for selection of filling
    - Option to show or hide waitbar added
%}
v6 : Suraj R Pawar, 6-12-2020
    - Added Qafilter structure as input
    - Structure houses the filtered signal, upper and lower thresholds for
    determining the stage of the cardiac cycle.
%}
    
    %% Argument handling
        if nargin < 14
            waitflag = 0;   % Don't display waitbar if the argument isn't passed in
        end
        
    %% UKF Parameters
        alpha = ukf_params(1);   
        kappa = ukf_params(2);
        beta = ukf_params(3);
         
    %% CVS Parameters
        [getparams, ~] = func_handle_parameters(parameters, version);
        
    %% Sweep Qa signal and index estimation windows
        steps = length(t);  % Total number of samples in measurement signals
        Plv = y(1,:);
        Pao = y(2,:);
        Qa = y(3,:); 
        Pr = getparams.Pr;
                
        tc = getparams.tc;  % Used for calculating normalized time (s)        
        j = 1;              % To index and store measurements in selected windows        
        Qcheck = 0;         % Should we start looking for Qa to drop below threshold ?        
        dt = t(2) - t(1);   
        
        % Uncomment for figures
        %{
        figure;
        plvplot = subplot(2,1,1);
        paoplot = subplot(2,1,2);
        qaplot = subplot(3,1,3);
        ax = [];
        
        axes(plvplot);
        hold on;
        plot(t,y(1,:),'k');
        title('Plv (mmHg)');
        ax = [ax, gca];
        
        axes(paoplot);
        hold on;
        plot(t,y(2,:),'k');
        title('Pao (mmHg)');
        ax = [ax, gca];
        
        axes(qaplot);
        hold on;
        plot(t,y(3,:),'k');
        title('Qa (mL/s)');
        %}
        
        % Selection of windows : This was with Plv based selection of
        % windows. Discarded this approach
        %{
        for i = 1:steps             
            tr = mod(t(i),tc);
            if tr >= 0 && tr <= dt
                rwt = 1;          
                marking = 0;
            end
            
            if rwt == 1
                if Plvfilt(i) > 100; % After R wave, wait for Plv to rise above 100 mmHg before setting Qcheck flag
                    Qcheck = 1;   
                    rwt = 0;
                    rst = 1;
                end                
            else
                if Qcheck == 1
                    if Plvfilt(i) <= 50                                    
                        marking = 1;
                        Qcheck = 0;
                    end                    
                end
            end
            
            if marking == 1                                  
                % filling window detected  
                Qa_selection(j) = Qa(i);
                Plv_selection(j) = Plv(i);
                Pao_selection(j) = Pao(i);
                Qvad_selection(j) = Qvad(i);
                tselect(j) = t(i);                      
                reset(j) = rst;             % 'reset' used to reset x and p for states during UKF
                j = j + 1;            
                if rst == 1
                    rst = 0;                                    
                end
            end                        
        end
        %}
        % Extract Qa filtered signal and thresholds
        Qafilt = Qafiltstruct.signal;
        upper_threshold = Qafiltstruct.upper_threshold;
        lower_threshold = Qafiltstruct.lower_threshold;
        
        for i = 1:steps 
            tr = mod(t(i),tc);
            if tr >= 0 && tr <= dt
                rwt = 1;                    
                marking = 0;
            end
            
            if rwt == 1
                if Qafilt(i) > upper_threshold % After R wave, wait for Qa to rise above 500 mL before setting Qcheck flag
                    Qcheck = 1;
                    rwt = 0;
                end                
            else
                if Qcheck == 1
                    if Qafilt(i) < lower_threshold
                        marking = 1;
                        Qcheck = 0;
                        rst = 1;                    
                    end
                end
            end
            
            if marking == 1                      
                % filling window detected  
                Qa_selection(j) = Qa(i);
                Plv_selection(j) = Plv(i);
                Pao_selection(j) = Pao(i);
                Qvad_selection(j) = Qvad(i);
                tselect(j) = t(i);      
                reset(j) = rst;             % 'reset' used to reset x and p for states during UKF
                j = j + 1;        
                if rst == 1
                    rst = 0;
                end
            end                        
        end
        
        % Uncomment for figure (match what was done with figure before loop)
        %{
        axes(plvplot);
        plot(tselect,Plv_selection,'r.');
        plot(tselect, reset.*80);        
        hold off;        
        
        axes(paoplot);
        plot(tselect,Pao_selection,'r.');
        hold off;
        
        linkaxes(ax,'x');
        
        axes(qaplot);
        plot(tselect,Qa_selection,'r.');
        hold off;
        %}            
        
    % Re-organize variables
        t = tselect;
        Qvad = Qvad_selection;        
        y = [Plv_selection; Pao_selection; Qa_selection];
        
    %% Variables and Initial conditions    
        % Count everything            
            dt = t(2) - t(1);                           % Time vector : [t0 : dt : tf]
            num_states = length(x0) + length(theta0);   % Number of states and parameters to estimate            
            num_process_noise_terms = size(q,1);
            num_meas_noise_terms = size(r,1);
            num_aug_states = num_states + num_process_noise_terms...
                             + num_meas_noise_terms;
            num_meas = size(y,1); 
            counts = [num_states; num_meas;...
                      num_process_noise_terms; num_meas_noise_terms]; % Information to pass on to Two Elem function file
            
        % Declare variables and set initial conditions        
            % Augmented state and covariance matrix
            steps_selection = numel(tselect);
            x_aug = zeros(num_aug_states, steps_selection);            
            x_aug0 = [x0; theta0; zeros(num_process_noise_terms,1);...
                          zeros(num_meas_noise_terms,1)];
            x_aug(:,1) = x_aug0;
            p_aug = zeros(num_aug_states, num_aug_states, steps_selection);
            QdB = diag(q)*dt;                                         % Noise covariance of brownian motion
            p_aug0 = blkdiag(p0, diag(QdB), r);
            p_aug(:,:,1) = p_aug0;
            
            % Estimated measurements
            yhat = zeros(num_meas,steps_selection);
            
            % Sigma points
            L = num_aug_states;
            lambda = alpha^2*(L + kappa) - L;
            wm_0 = lambda/(L + lambda);        
            wc_0 = (lambda/(L + lambda)) + (1 - alpha^2 + beta); 
            wm_i = 1/(2*(lambda + L));          
            wc_i = wm_i;                      
            wm_all = [wm_0, repmat(wm_i,1,2*L)];
            wc_all = [wc_0, repmat(wc_i,1,2*L)];
            num_sigmas = 2*L + 1;
            sigmas = zeros(num_aug_states, num_sigmas, steps_selection);                                                  
            
            % Frequency of console display            
            if waitflag == 1 
                f = waitbar(0,'UKF'); 
                console_freq = floor(steps_selection/10);
            end
            
    %% UKF   
        for i = 2:steps_selection
            % Console out
                if waitflag == 1
                    if mod(i,console_freq) == 0                    
                        waitbar((i/steps_selection),f,'UKF');
                    end
                end
            
            % Reset UKF at beginning of ejection window
                if reset(i-1) == 1
                    x_aug([1:2],i-1) = x_aug0([1:2]);             
                    p_aug([1:2],[1:2],i-1) = p_aug0([1:2],[1:2]);
                end
                
            % Extract previous augmented state and cov matrix
                x_previous = x_aug(:,i-1);
                p_previous = squeeze(p_aug(:,:,i-1));
                
            % Prepare sigma points
                sqrt_mat = chol((L + lambda)*p_previous).';
                sigmas(:,1,i-1) = x_previous;
                sigmas(:,[2 : L+1], i-1) = x_previous + sqrt_mat;
                sigmas(:,[L+2 : end], i-1) = x_previous - sqrt_mat;                
                
            % Time Update - propagate sigma points    
                %{
                  Its ok to use the state equations for ejection/filling, because
                  Qa/Qm will appropriately set to 0 in these state equations
                %}
                [sigmas_x_new, est_meas] = func_TwoElem_Filling(t(i-1), dt, squeeze(sigmas(:,:,i-1)),...
                                             Qvad(i-1), parameters, counts, version);
                
            % Prior mean and cov of sigma points
                mean_prior = sum((wm_all.*sigmas_x_new),2);                    
                for j = 1 : (2*L + 1)
                    if j == 1
                        cov_prior = zeros(size(p0));
                    end
                    cov_prior = cov_prior  + wc_all(j)*(sigmas_x_new(:,j) - mean_prior)*(sigmas_x_new(:,j)-mean_prior).';                        
                end                    
                                                                            
            % Measurement update
                mean_prior_meas = sum((wm_all.*est_meas),2);
                for kk = 1 : (2*L + 1)
                    if kk == 1
                        cov_meas_prior = zeros(num_meas);
                        cov_cross_prior = zeros(num_states,num_meas);
                    end
                    cov_meas_prior = cov_meas_prior + wc_all(kk)*(est_meas(:,kk) - mean_prior_meas)*(est_meas(:,kk) - mean_prior_meas).';
                    cov_cross_prior = cov_cross_prior + wc_all(kk)*(sigmas_x_new(:,kk) - mean_prior)*(est_meas(:,kk) - mean_prior_meas).';
                end                
                ukf_gain = cov_cross_prior*inv(cov_meas_prior);
                
                mean_post = mean_prior + ukf_gain*(y(:,i) - mean_prior_meas);
                cov_post = cov_prior - ukf_gain*cov_meas_prior*ukf_gain.';                       
                
            % Update the augmented state vector
                x_aug([1:num_states],i) = mean_post;
                p_aug(:,:,i) = blkdiag(cov_post, diag(QdB), r);   
                
            % Update measurements
                yhat(:,i) = func_update_measurements(t(i), mean_post, mean_prior_meas, parameters, 3, version);                                                                  
        end
        
        if waitflag == 1
            close(f);
        end
        
    %% Generate Outputs
        xhat = x_aug([1:num_states],:);        
end
