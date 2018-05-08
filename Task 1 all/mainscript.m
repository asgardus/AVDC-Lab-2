%% Clearing workspace
clear all
clc

%% Data Loading

global WA_VBOX vx_VBOX Time_data file

% file_list = ["crc","sla","swd","step"];
% for b=1:length(file_list)
%     file = file_list(b);
file = 'crc'; %specify file here according to below specified legend
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

cf_mat=80e3;      %input range or single value
cr_mat=85e3;      %input range or single value
T_w_mat=1:0.5:10;        %input range or single value:
for i=1:length(cf_mat)
        Cf=cf_mat(i);
        for j=1:length(cr_mat)
            Cr=cr_mat(j);            
            for k=1:length(T_w_mat)
                T_w=T_w_mat(k);
                sim('Estimator_Model');
                % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
                [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_washout_main.data, Beta_VBOX);
                mse(k)=e_beta_mean;
                max(k)=e_beta_max;
            end
        end
end
figure(2);
title('Filter coefficient tuning results')
subplot(1,2,1)
plot(T_w_mat,mse,'DisplayName',file)
xlabel('Filter coefficient T')
ylabel('Mean Squared Error (rad)')
ylim([0 0.004])
grid on
hold on
legend('show')
subplot(1,2,2);
plot(T_w_mat,max,'DisplayName',file)
xlabel('Filter coefficient T')
ylabel('Max Error (rad)')
grid on
hold on
legend('show')

%% Task 1-d

% cf_mat=80e3;      %input range or single value
% cr_mat=85e3;      %input range or single value
% T_w_mat=3.5;        %input range or single value:
% for i=1:length(cf_mat)
%         Cf=cf_mat(i);
%         for j=1:length(cr_mat)
%             Cr=cr_mat(j);            
%             for k=1:length(T_w_mat)
%                 T_w=T_w_mat(k);
%                 sim('Estimator_Model');
%                 % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
%                 [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_washout_main.data, Beta_VBOX);
%                 mse(1)=e_beta_mean;
%                 max(1)=e_beta_max;
%                 % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
%                 [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_bicycle_sim_main.data, Beta_VBOX);
%                 mse(2)=e_beta_mean;
%                 max(2)=e_beta_max;
%                 % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
%                 [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_integrator_main.data, Beta_VBOX);
%                 mse(3)=e_beta_mean;
%                 max(3)=e_beta_max;
%             end
%         end
% end
% figure(2);
% x=1:3;
% subplot(1,2,1)
% title('Performance Evaluation-Mean Squared Error')
% plot(x,mse,'DisplayName',file)
% ylabel('Mean Squared Error (rad)')
% grid on
% hold on
% legend('show')
% 
% subplot(1,2,2);
% title('Performance Evaluation-Max Error')
% plot(x,max,'DisplayName',file)
% ylabel('Max Error (rad)')
% grid on
% hold on
% legend('show')
% end
%% Task 1-g
% m_mat=9.5;      %input range or single value
% c_mat=1;        %input range or single value
% Tau = 0.3;
% T_w = 2;
% count=1;
% for i=1:length(m_mat)
%         m=m_mat(i);
%         for j=1:length(c_mat)
%             c=c_mat(j);
%             sim('Estimator_Model_T_variable');
%             % CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE             
%             [e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(beta_washout_main.data, Beta_VBOX);
%             mse(count)=e_beta_mean;
%             max(count)=e_beta_max;         
%             counter(count)=count;
%             count=count+1;
%         end
% end
% counter(1)=1;
% subplot(1,2,1)
% title('Performance Evaluation-Mean Squared Error')
% plot(counter,mse,'DisplayName',file)
% ylabel('Mean Squared Error (rad)')
% grid on
% hold on
% legend('show')
% 
% subplot(1,2,2);
% title('Performance Evaluation-Max Error')
% plot(counter,max,'DisplayName',file)
% ylabel('Max Error (rad)')
% grid on
% hold on
% legend('show')
% 
% mse_error(:,b)=mse;
% max_error(:,b)=max;
% % end