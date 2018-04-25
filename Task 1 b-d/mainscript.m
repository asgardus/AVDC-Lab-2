Init_for_washout_filter;

global WA_VBOX vx_VBOX Time

WA_VBOX = SWA_VBOX./Ks;
WA_VBOX_mat = [Time WA_VBOX];
yawRate_VBOX_mat = [Time yawRate_VBOX];
vx_VBOX_mat = [Time vx_VBOX];
ay_VBOX_mat = [Time ay_VBOX];
Beta_VBOX_mat = [Time Beta_VBOX];
vx = vx_VBOX;
t = Time;
x0 = [-0.0103 0.1244 0.00052 0];
T = Time(end);
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
Cf_tune = 10000:10000:110000;
Cr_tune = 10000:10000:110000;

for i=1:length(Cf_tune)
    Cf = Cf_tune(i);
    for j=1:length(Cr_tune)
        Cr = Cr_tune(j);
        sim('Estimator_Model');
    end
end

%% Plot results
figure(1);
plot(Time, Beta_VBOX, Time, beta_bicycle_sim_main.Data);