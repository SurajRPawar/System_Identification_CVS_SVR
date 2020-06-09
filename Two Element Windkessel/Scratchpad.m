clear all; close all; clc;

num_aug_states = 4;
num_meas_ejection = 3;

Lf = sym(zeros(num_aug_states*num_meas_ejection,1));

for i = 1:num_aug_states
    if i == 1
        start_index = 1;
    else
        start_index = (i-1)*(num_meas_ejection) + 1;
    end
    stop_index = i*num_meas_ejection;
    Lf(start_index : stop_index) = i*ones(num_meas_ejection,1);
end

Lf