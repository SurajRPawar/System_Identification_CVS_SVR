function [getparams, paramdots] = func_handle_parameters(parameters, version, sigmas, x)
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
v1 : Suraj R Pawar, 5-28-2020
    - Initialize
%}
    
    num_states = 2; % Number of states in the two element windkessel model (2)
    switch version
        case 1
            % Estimation of B and Emax
            getparams.Rsvr = parameters(1);
            getparams.Cs = parameters(2);
            getparams.Pr = parameters(3);
            getparams.Rv = parameters(4);
            getparams.HR = parameters(6);
            getparams.tvc = parameters(7);
            getparams.tc = parameters(8);
            getparams.A = parameters(9);
            if nargin == 3
                % Sigma points passed
                getparams.B = sigmas((num_states+1),:);
                getparams.E = sigmas((num_states+2),:);        
                paramdots = zeros(2,size(sigmas,2));
            elseif nargin == 4
                % Last best estimate of the augmented state passed
                getparams.B = x(3);
                getparams.E = x(4);
                paramdots = zeros(2,1);
            elseif nargin == 2
                % UKF has not yet begun
                getparams.B = 0;
                getparams.E = 0;
                paramdots = zeros(2,1);
            end
        case 2             
            % Estimation of A
            getparams.Rsvr = parameters(1);
            getparams.Cs = parameters(2);
            getparams.Pr = parameters(3);
            getparams.Rv = parameters(4);
            getparams.HR = parameters(6);
            getparams.tvc = parameters(7);
            getparams.tc = parameters(8);
            getparams.B = parameters(10);
            getparams.E = parameters(11);
            if nargin == 3
                % Sigma points passed
                getparams.A = sigmas((num_states+1),:);                
                paramdots = zeros(1,size(sigmas,2));
            elseif nargin == 4
                % Last best estimate of the augmented state passed
                getparams.A = x(3);                
                paramdots = zeros(1,1);
            elseif nargin == 2
                % UKF has not yet begun
                getparams.A = 0;                
                paramdots = zeros(1,1);
            end              
        case 3
            % Estimation of B and Emax
            getparams.Rsvr = parameters(1);
            getparams.Cs = parameters(2);
            getparams.Pr = parameters(3);
            getparams.Rv = parameters(4);
            getparams.HR = parameters(6);
            getparams.tvc = parameters(7);
            getparams.tc = parameters(8);            
            if nargin == 3
                % Sigma points passed
                getparams.A = sigmas((num_states+1),:);
                getparams.B = sigmas((num_states+2),:);
                getparams.E = sigmas((num_states+3),:);        
                paramdots = zeros(3,size(sigmas,2));
            elseif nargin == 4
                % Last best estimate of the augmented state passed
                getparams.A = x(3);
                getparams.B = x(4);
                getparams.E = x(5);
                paramdots = zeros(3,1);
            elseif nargin == 2
                % UKF has not yet begun
                getparams.A = 0;
                getparams.B = 0;
                getparams.E = 0;
                paramdots = zeros(3,1);
            end            
    end
       
end