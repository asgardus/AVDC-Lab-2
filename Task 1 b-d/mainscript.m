Init_for_washout_filter;

global WA_VBOX vx_VBOX Time

WA_VBOX = SWA_VBOX./Ks
SWA_VBOX_mat = [Time WA_VBOX];
yawRate_VBOX_mat = [Time yawRate_VBOX];
vx_VBOX_mat = [Time vx_VBOX];
ay_VBOX_mat = [Time ay_VBOX];
Beta_VBOX_mat = [Time Beta_VBOX];
vx = vx_VBOX;
t = Time;
x0 = [-0.0103 0.1244 0.00052 0];
T = Time(end)+0.01;
T_w = 1;

%% MATLAB Bicycle Model Estimator
% [timeout,xout]=ode45('bicycle_estimator_mat', t, x0);
% vy = xout(:,1);
% beta_bicycle_mat = [Time atan(vy./vx)];
% 
% %% Simulink Integrator Estimator
% sim('integrator_estimator');
% 
% %% Simulink Bicycle Model Estimator
% sim('bicycle_estimator');
% 
% %% Simulink Washout Filter Estimator
% sim('washout_estimator');

%% Simulink Full Estimator Model
sim('Estimator_Model');

%% Plot results
len = [length(beta_integrator_main.Data) length(beta_bicycle_sim_main.Data) length(beta_washout_main.Data) length(Beta_VBOX)];
figure(1);
plot(Time, Beta_VBOX, Time, beta_integrator_main.Data, Time, beta_bicycle_sim_main.Data, Time, beta_washout_main.Data);