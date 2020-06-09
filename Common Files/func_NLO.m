%--------------------------------------------------------------------------
% Function to calculate the nonlinear observability rank. All inputs must
% be symbolics ! This function uses the formulation due to Hermann and
% Krener.
% Inputs = f : Vector field as function of t, x, u
%          y : Meaurement vector
%          x : State vector
% Outputs = r : Rank of nonlinear observability matrix
%           dLf : Annhilator of the observability matrix
%
% Suraj R Pawar, 4-17-2020
%--------------------------------------------------------------------------

function [r, dLf] = func_NLO(x,f,y)

    % Define variables
        num_states = length(f);
        num_meas = length(y);
        num_elems_lie = num_states*num_meas; % No. of columns x No. of rows = Total number of elements in the column vector formed by stacking all repeated Lie Derivatives
        Lf = sym(zeros(num_elems_lie,1));   % Variable to store all repeated Lie derivatives. This vector will be used to form the co-observability space
        dLf = sym(zeros(num_elems_lie,num_states));
        
    % Calculate Observability space, and annhilator of Obs space
        for i = 1:num_states
            if i == 1
                start_index = 1;
                stop_index = i*num_meas;
            else
                start_index = (i-1)*num_meas + 1;
                start_index_previous = (i-2)*num_meas + 1;
                stop_index = i*num_meas;
                stop_index_previous = (i-1)*num_meas;
            end
            
            
            if i == 1
                    Lf(start_index:stop_index) = y;
                else
                    Lf(start_index:stop_index) = jacobian(Lf(start_index_previous:stop_index_previous),x)*f;
            end
            
            dLf([start_index:stop_index],:) = jacobian(Lf(start_index:stop_index),x);
        end
        
    % Calculate rank of annhilator of observation space. 
    r = rank(dLf);    
end