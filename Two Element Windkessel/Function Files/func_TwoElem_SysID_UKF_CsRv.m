function [xhat, yhat, p_aug] = func_TwoElem_SysID_UKF_CsRv(t, y, u, x0, theta0,...
                                                           p0, q, r, parameters,...
                                                           ukf_params, version, waitflag)        
% UKF for Cs and Rv estimation using Two Element Windkessel model
%{
----------------------------- DESCRIPTION ---------------------------------
Take the interpolated measurements of Pao, and apply UKF to the single
order ODE for estimation of Cs and Pr

----------------------------- ASSUMPTIONS ---------------------------------
- We know SVR from the equation (MAP - CVP)/CO
- We know A, B, Emax from a separate UKF
- We know Rv from mean dP / mean Qa

----------------------------- INPUTS --------------------------------------
t          : Time vector in the format [t0 : dt: tf]
y          : Measurement signals of Pao
Qao        : Aortic flow rate input (QVAD + Qa) (mL/s)
x0         : Initial guess for states
theta0     : initial guess for parameters [Cs; Pr]
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
waitflag (o) : Wait flag (0 or 1) to display GUI waitbar. 1 = display
               waitbar

----------------------------- OUTPUTS -------------------------------------
xhat    : Estimated states
theta   : Estimated parameters
p_aug   : Covariances for the augmented state [x; theta; v; n]

----------------------------- VERSION -------------------------------------
v1 : Suraj R Pawar, 6-11-2020
    - Initialize
%}
    
    %% Argument handling
        if nargin < 12
            waitflag = 0;   % Don't display waitbar if the argument isn't passed in
        end
        
    %% UKF Parameters
        alpha = ukf_params(1);   
        kappa = ukf_params(2);
        beta = ukf_params(3);             
        
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
            steps = numel(t);
            x_aug = zeros(num_aug_states, steps);            
            x_aug0 = [x0; theta0; zeros(num_process_noise_terms,1);...
                          zeros(num_meas_noise_terms,1)];
            x_aug(:,1) = x_aug0;
            p_aug = zeros(num_aug_states, num_aug_states, steps);
            QdB = diag(q)*dt;                                         % Noise covariance of brownian motion
            p_aug0 = blkdiag(p0, diag(QdB), r);
            p_aug(:,:,1) = p_aug0;
            
            % Estimated measurements
            yhat = zeros(num_meas,steps);
            
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
            sigmas = zeros(num_aug_states, num_sigmas, steps);                                                  
            
            % Frequency of console display            
            if waitflag == 1 
                f = waitbar(0,'UKF'); 
                console_freq = floor(steps/10);
            end
            
    %% UKF   
        for i = 2:steps
            % Console out
                if waitflag == 1
                    if mod(i,console_freq) == 0                    
                        waitbar((i/steps),f,'UKF');
                    end
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
                [sigmas_x_new, est_meas] = func_TwoElem_CsRv(t(i-1), dt, squeeze(sigmas(:,:,i-1)),...
                                             u(:,i-1), parameters, counts, version);
                
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
                yhat(:,i) = mean_prior_meas;
        end
        
        if waitflag == 1
            close(f);
        end
        
    %% Generate Outputs
        xhat = x_aug([1:num_states],:);
end
