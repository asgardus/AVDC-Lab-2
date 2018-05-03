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
% Cf_tune = 10000:10000:110000;
% Cr_tune = 10000:10000:110000;
% 
% for i=1:length(Cf_tune)
%     Cf = Cf_tune(i);
%     for j=1:length(Cr_tune)
%         Cr = Cr_tune(j);
%         sim('Estimator_Model');
%     end
% end

%% Mean squared error
% for i=1:100
%     j(i)=i;
%         tw(i)=i*0.005;
%         T_w=i*0.05;
%          sim('Estimator_Model');
% %         err(i) = immse(beta_washout_main.data,Beta_VBOX);
%        % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE
% %--------------------------------------------------------- 
% [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_washout_main.data, Beta_VBOX);
% disp(' ');
% % fprintf('The MSE of Beta estimation is: %d \n',e_beta_mean(i));
% % fprintf('The Max error of Beta estimation is: %d \n',e_beta_max);
% mse(i)=e_beta_mean
% max(i)=e_beta_max
% end
% figure (2)        
% plot(tw,mse)
% % plot(tw,e_beta_max)
% hold on
% grid on

%% Yaw acceleration dependent Time function
m = 3;
c = 0.02;
Tau = 0.3;
sim('Estimator_Model_T_variable');

%% Plot results
figure(1);
plot(Time, Beta_VBOX, Time, beta_bicycle_sim_main.Data);