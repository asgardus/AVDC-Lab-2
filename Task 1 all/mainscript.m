%% Clearing workspace
clear all
clc

%% Data Loading

global WA_VBOX vx_VBOX Time_data file

file = 'sla'; %specify file here according to below specified legend
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

%% Task 1-a, 1-b and 1-c

cf_mat=80e3;      %select range or single value
cr_mat=85e3;      %select range or single value
T_w_mat=3.5;            %select range or single value
for i=1:length(cf_mat)
        Cf=cf_mat(i);
        for j=1:length(cr_mat)
            Cr=cr_mat(j);            
            for k=1:length(T_w_mat)
                T_w=T_w_mat(k);
                sim('Estimator_Model');
                % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
                [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_washout_main.data, Beta_VBOX);
                mse(i)=e_beta_mean;
                max(i)=e_beta_max;
            end
        end
end
figure(2);        
plot(mse)
grid on
figure(3);
plot(max)
grid on

%% Task 1-g
% m_mat=0.1:0.5:10;
% c_mat=eps:0.1;
% Tau = 0.3;
% T_w = 3.5;
% for i=1:length(m_mat)
%         m=m_mat(i);
%         for j=1:length(c_mat)
%             c=c_mat(j);
%             sim('Estimator_Model_T_variable');
%             % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
%             [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_washout_main.data, Beta_VBOX);
%             mse(i)=e_beta_mean;
%             max(i)=e_beta_max;
%         end
% end
% figure(2);        
% plot(mse)
% grid on
% figure(3);
% plot(max)
% grid on

%% Plot results
% figure(4);
% plot(Time, Beta_VBOX, Time, beta_bicycle_sim_main.Data);