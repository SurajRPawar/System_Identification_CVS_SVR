function [xhat, yhat, p_aug, tselect] = func_TwoElem_SysID_UKF_BE_Ejection(t, y, Qvad, x0, theta0,...
                                                                              p0, q, r, parameters,...
                                                                              ukf_params, version, Qafiltstruct, waitflag)        
% UKF for B and Emax estimation using Two Element Windkessel model
%{
----------------------------- DESCRIPTION ---------------------------------
Take measurement signals, detect isovolumic contraction and ejection and 
apply UKF only to these windows. At the beginning of each window, we reset
the states to their initial value, and reset the covariances corresponsing
to the states to the initial covariance. By states, we refer to [Vbar; Ps].
The estimates and covariances of the parameters being estimated will carry
forward.

----------------------------- ASSUMPTIONS ---------------------------------
- We know SVR from the equation (MAP - CVP)/CO
- We know Cs, Pr, A, Rv
- Only unknown parameters to estimate are B, Emax

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
Qafildt    : Filtered aortic flow signal, used for selection of filling
             stage
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
    - Updated comments
    - Removed input 'stage' from function definition because it was not
    being used anywhere
v5 : Suraj R Pawar, 6-10-2020
    - Added waitbar option
    - Cleaned up comments and description
v6 : Suraj R Pawar, 6-10-2020
    - Updated Qa threshold to stop marking ejection.
    - New Threshold is 0
v7 : Suraj R Pawar, 6-12-2020
    - Added Qafilter structure as input
    - Structure houses the filtered signal, upper and lower thresholds for
    determining the stage of the cardiac cycle.
%}
v8 : Suraj R Pawar, 6-14-2020
    - Fixed number of arguments being checked to set waitflag
%}
    
    %% Argument handling
        if nargin < 13
            waitflag = 0;   % If no input provided for waitflag, set it to 0
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
        
        marking = 0;        % Indicates when marking is in progress   
        tc = getparams.tc;  % Used for calculating normalized time (s)
        j = 1;              % To index and store measurements in selected windows
        rwt = 0;            % R wave trigger flag
        Qcheck = 0;         % Should we start looking for Qa to drop below threshold ?
        dt = t(2) - t(1);   
        
        % Uncomment for figures
        %{
        figure;
        plvplot = subplot(3,1,1);
        paoplot = subplot(3,1,2);
        qaplot = subplot(3,1,3);
        
        axes(plvplot);
        hold on;
        plot(t,y(1,:),'k');
        title('Plv (mmHg)');
        
        axes(paoplot);
        hold on;
        plot(t,y(2,:),'k');
        title('Pao (mmHg)');
        
        axes(qaplot);
        hold on;
        plot(t,y(3,:),'k');
        title('Qa (mL/s)');
        %}
        
        % Extract the right stage from measurements
        
        % Extract Qa filtered signal and thresholds
        Qafilt = Qafiltstruct.signal;
        upper_threshold = Qafiltstruct.upper_threshold;
        lower_threshold = Qafiltstruct.lower_threshold;
        
        for i = 1:steps 
            tr = mod(t(i),tc);
            if tr >= 0 && tr <= dt
                rwt = 1;    
                rst = 1;
            else
                rst = 0;
            end
            
            if rwt == 1
                if Qafilt(i) > upper_threshold % After R wave, wait for Qa to rise above 500 mL before setting Qcheck flag
                    Qcheck = 1;
                    rwt = 0;
                end
                marking = 1;
            else
                if Qcheck == 1
                    if Qafilt(i) > lower_threshold
                        marking = 1;
                    else
                        marking = 0;
                        Qcheck = 0;
                    end
                end
            end
            
            if marking == 1                      
                % Ejection window detected  
                Qa_ejection(j) = Qa(i);
                Plv_ejection(j) = Plv(i);
                Pao_ejection(j) = Pao(i);
                Qvad_ejection(j) = Qvad(i);
                tselect(j) = t(i);      
                reset(j) = rst;             % 'reset' used to reset x and p for states during UKF
                j = j + 1;            
            end                        
        end
        
        % Uncomment for figure (match what was done with figure before loop)
        %{
        axes(plvplot);
        plot(tselect,Plv_ejection,'r.');
        plot(tselect, reset.*80);
        hold off;        
        
        axes(paoplot);
        plot(tselect,Pao_ejection,'r.');
        hold off;
        
        axes(qaplot);
        plot(tselect,Qa_ejection,'r.');
        hold off;
        %}        
%           1;
        % Re-organize variables
        t = tselect;
        Qvad = Qvad_ejection;
        y = [Plv_ejection; Pao_ejection; Qa_ejection];
        
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
            steps_ejection = numel(tselect);
            x_aug = zeros(num_aug_states, steps_ejection);            
            x_aug0 = [x0; theta0; zeros(num_process_noise_terms,1);...
                          zeros(num_meas_noise_terms,1)];
            x_aug(:,1) = x_aug0;
            p_aug = zeros(num_aug_states, num_aug_states, steps_ejection);
            QdB = diag(q)*dt;                                         % Noise covariance of brownian motion
            p_aug0 = blkdiag(p0, diag(QdB), r);
            p_aug(:,:,1) = p_aug0;
            
            % Estimated measurements
            yhat = zeros(num_meas,steps_ejection);
            
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
            sigmas = zeros(num_aug_states, num_sigmas, steps_ejection);                                                  
            
    %% UKF   
        % Frequency of console display            
        if waitflag == 1 
            f = waitbar(0,'UKF'); 
            console_freq = floor(steps_ejection/10);
        end
        
        for i = 2:steps_ejection
            % Console out
                if waitflag == 1
                    if mod(i,console_freq) == 0                    
                        waitbar((i/steps_ejection),f,'UKF');
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
                  Its ok to use the state equations for ejection, because
                  Qa will appropriately set to 0 in these state equations
                %}
                [sigmas_x_new, est_meas] = func_TwoElem_Ejection(t(i-1), dt, squeeze(sigmas(:,:,i-1)),...
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
                yhat(:,i) = func_update_measurements(t(i), mean_post, mean_prior_meas, parameters, 2, version);                                                                     
        end
 
        if waitflag == 1
            close(f);
        end
        
    %% Generate Outputs
        xhat = x_aug([1:num_states],:);        
end
