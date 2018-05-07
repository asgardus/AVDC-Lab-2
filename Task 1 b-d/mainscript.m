%% Clearing workspace
clear all
clc

%% Data Loading

global WA_VBOX vx_VBOX Time_data file

file = 'stand'; %specify file here according to below specified legend
%'stand' - standstill       'crc' - circle test left        'swd' - sine dwell
%'sla' - slaloms            'step' - step steer

Init_for_washout_filter;

switch file
    case 'sla'
        Time1 = Time_data(1):0.01:Time_data(end)+0.01;
    case 'stand'
        Time1 = Time_data(1):0.01:Time_data(end)+0.01;
    otherwise
        Time1 = Time_data(1):0.01:Time_data(end);
end
Time = Time1';
WA_VBOX = SWA_VBOX./Ks;
WA_VBOX_mat = [Time WA_VBOX];
yawRate_VBOX_mat = [Time yawRate_VBOX];
vx_VBOX_mat = [Time vx_VBOX];
ay_VBOX_mat = [Time ay_VBOX];
Beta_VBOX_mat = [Time Beta_VBOX];
roll_angle_VBOX_mat = [Time roll_angle_VBOX];
vx = vx_VBOX;
t = Time;
x0 = [-0.0103 0.1244 0.00052 0];
T = Time(end)-Time(1);
T_w = 3.5;

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
m_mat=0.1:0.5:10;
c_mat=eps:0.1;
Tau = 0.3;
for i=1:length(m_mat)
        m=m_mat(i);
        for j=1:length(c_mat)
            c=c_mat(j);
            sim('Estimator_Model_T_variable');
            % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
            [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_washout_main.data, Beta_VBOX);
            mse(i)=e_beta_mean;
            max(i)=e_beta_max;
        end
end
figure(2);        
plot(mse)
grid on
figure(3);
plot(max)
grid on

%% Yaw acceleration dependent Time function
% m = 3;
% c = 0.02;
% sim('Estimator_Model_T_variable');

%% Plot results
% figure(1);
% plot(Time, Beta_VBOX, Time, beta_bicycle_sim_main.Data);