%% Perform checks on the estimated parameters
% if any(sigmas_x_new(3,:) < 0)
%     sigmas_x_new(3,find(sigmas_x_new(3,:) < 0)) = 0;
% end
% if any(sigmas_x_new(4,:) < 0)
%     sigmas_x_new(4,find(sigmas_x_new(4,:) < 0)) = 0;
% end
% if any(sigmas_x_new(5,:) < 0)
%     sigmas_x_new(5,find(sigmas_x_new(5,:) < 0)) = 0;
% end