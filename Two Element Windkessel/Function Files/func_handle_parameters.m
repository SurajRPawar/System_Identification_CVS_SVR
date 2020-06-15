function [getparams, paramdots] = func_handle_parameters(parameters, version, counts, sigmas, x)
%% Extract parameters and set parameter state equations based on versions
%{
----------------------------- Description ---------------------------------
Based on the version of the experiment being run, extract the true and
estimated values of the parameters. Further, set the dynamic equations of
parameters being estimated accordingly. For a description of the versions,
open the MATLAB script file 'experiment_versions.m'

----------------------------- Inputs --------------------------------------
parameters : Vector containing true value of all parameters
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
version    : Version number of the experiment
counts     : 1. num_states : Number of states and parameters to be
                estimated
             2. num_meas : Number of measurements
             3. num_process_noise_terms : Number of process noise terms
             4. num_meas_noise_terms : Number of meas noise terms
sigmas (o) : Sigma points (optional)
             [sigma_x; -- States and parameters to be estimated
              sigma_v; -- Process noise terms
              sigma_n] -- Measurement noise terms
x (o)      : Estimted states (optional)

----------------------------- Outputs -------------------------------------
getparams  : Parameters to be used in model file that contains state
             equations
paramdots  : Dynamics of parameters to be estimated

------------------------------- Versions ----------------------------------
%{
v1 : Suraj R Pawar, 5-28-2020
    - Initialize
v2 : Suraj R Pawar, 6-10-2020
    - Clean up comments describing new cases
v3 : Suraj R Pawar, 6-11-2020
    - Add counts to the inputs so num_states doesn't have to be hardcoded
%}
v4 : Suraj R Pawar, 6-14-2020
    - Add case which handles only Cs estimation
    - Add variable of num_params that decides number of rows in paramdots
%}
        
    switch version
        case 1
            % Estimation of B and Emax
            num_params = 2;  % Number of parameters being estimated in this case            
            getparams.Rsvr = parameters(1);
            getparams.Cs = parameters(2);
            getparams.Pr = parameters(3);
            getparams.Rv = parameters(4);
            getparams.HR = parameters(6);
            getparams.tvc = parameters(7);
            getparams.tc = parameters(8);
            getparams.A = parameters(9);
            if nargin == 4
                % Sigma points passed
                num_states = counts(1);                
                getparams.B = sigmas((num_states-1),:);
                getparams.E = sigmas((num_states),:);        
                paramdots = zeros(num_params,size(sigmas,2));
            elseif nargin == 5
                % Last best estimate of the augmented state passed
                getparams.B = x(3);
                getparams.E = x(4);
                paramdots = zeros(num_params,1);
            elseif nargin == 2
                % UKF has not yet begun
                getparams.B = 0;
                getparams.E = 0;
                paramdots = zeros(num_params,1);
            end
        case 2
            % Estimation of A
            num_params = 1;
            getparams.Rsvr = parameters(1);
            getparams.Cs = parameters(2);
            getparams.Pr = parameters(3);
            getparams.Rv = parameters(4);
            getparams.HR = parameters(6);
            getparams.tvc = parameters(7);
            getparams.tc = parameters(8);
            getparams.B = parameters(10);
            getparams.E = parameters(11);
            if nargin == 4
                % Sigma points passed
                num_states = counts(1);
                getparams.A = sigmas((num_states),:);                
                paramdots = zeros(num_params,size(sigmas,2));
            elseif nargin == 5
                % Last best estimate of the augmented state passed
                getparams.A = x(3);                
                paramdots = zeros(num_params,1);
            elseif nargin == 2
                % UKF has not yet begun
                getparams.A = 0;                
                paramdots = zeros(num_params,1);
            end              
        case 3
            % Estimation of B, Emax and A
            num_params = 3;
            getparams.Rsvr = parameters(1);
            getparams.Cs = parameters(2);
            getparams.Pr = parameters(3);
            getparams.Rv = parameters(4);
            getparams.HR = parameters(6);
            getparams.tvc = parameters(7);
            getparams.tc = parameters(8);            
            if nargin == 4
                % Sigma points passed
                num_states = counts(1);
                getparams.A = sigmas((num_states-2),:);
                getparams.B = sigmas((num_states-1),:);
                getparams.E = sigmas((num_states),:);        
                paramdots = zeros(num_params,size(sigmas,2));
            elseif nargin == 5
                num_states = counts(1);
                % Last best estimate of the augmented state passed
                getparams.A = x(3);
                getparams.B = x(4);
                getparams.E = x(5);
                paramdots = zeros(num_params,1);
            elseif nargin == 2
                % UKF has not yet begun
                getparams.A = 0;
                getparams.B = 0;
                getparams.E = 0;
                paramdots = zeros(num_params,1);
            end    
        case 4
            % Estimation of Cs and Pr
            num_states = counts(1);
            num_params = 2;
            getparams.Rsvr = parameters(1);
            getparams.Cs = sigmas((num_states-1),:);
            getparams.Pr = sigmas((num_states),:);
            paramdots = zeros(num_params,size(sigmas,2));
        case 5
            % Estimation of Cs only
            num_states = counts(1);
            num_params = 1;
            getparams.Rsvr = parameters(1);
            getparams.Cs = sigmas((num_states),:);
            getparams.Pr = parameters(3);
            paramdots = zeros(num_params,size(sigmas,2));
        case 6
            % Estimation of A and Rv
            num_params = 2;
            getparams.Rsvr = parameters(1);
            getparams.Cs = parameters(2);
            getparams.Pr = parameters(3);            
            getparams.HR = parameters(6);
            getparams.tvc = parameters(7);
            getparams.tc = parameters(8);
            getparams.B = parameters(10);
            getparams.E = parameters(11);
            if nargin == 4
                % Sigma points passed
                num_states = counts(1);
                getparams.A = sigmas((num_states - 1),:);   
                getparams.Rv = sigmas((num_states),:);   
                paramdots = zeros(num_params,size(sigmas,2));
            elseif nargin == 5
                % Last best estimate of the augmented state passed
                getparams.A = x(3);     
                getparams.Rv = x(4);     
                paramdots = zeros(num_params,1);
            elseif nargin == 2
                % UKF has not yet begun
                getparams.A = 0;         
                getparams.Rv = 0;
                paramdots = zeros(num_params,1);
            end              
    end
       
end