%----------------------------------------------------------------
% Template created for the course SD2231 by Mikael Nybacka 2013
% Following file is the start file for the state estimation using
% Uncented Kalman Filter (UKF).
%----------------------------------------------------------------
clear all;
close all;
clc;
addpath('scripts')
addpath('logged_data')
disp(' ');
 
%% Data Loading

global Cf Cr mass Iz WA_VBOX vx_VBOX Time_data file dt delta Ratio SWA_VBOX

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
% WA_VBOX = SWA_VBOX./Ks;
% WA_VBOX_mat = [Time WA_VBOX];
% yawRate_VBOX_mat = [Time yawRate_VBOX];
% vx_VBOX_mat = [Time vx_VBOX];
% ay_VBOX_mat = [Time ay_VBOX];
% Beta_VBOX_mat = [Time Beta_VBOX];
% roll_angle_VBOX_mat = [Time roll_angle_VBOX];
% vx = vx_VBOX;
% t = Time;
% x0 = [-0.0103 0.1244 0.00052 0];
% T = Time(end)-Time(1);
%% UKF
% Set global variables so that they can be accessed from other matlab
% functions and files
% global lf lr Cf Cr mass Iz vbox_file_name dt SWA_VBOX Ratio delta

%----------------------------
% LOAD DATA FROM VBOX SYSTEM
%----------------------------
% vbox_file_name='logged_data/Lunda_test_140411/Stand_Still_no2.VBO'; %stand still logging, engine running
% vbox_file_name='logged_data/Lunda_test_140411/Circle_left_R13m_no2.VBO'; %circle test left, roughly 13m in radius
% vbox_file_name='logged_data/Lunda_test_140411/Slalom_35kph.VBO'; %slalom entry to the left @ first cone, 35kph
% vbox_file_name='logged_data/Lunda_test_140411/Step_Steer_left_80kph.VBO'; %Step steer to the left in 80kph
% vbox_file_name='logged_data/Lunda_test_140411/SWD_80kph.VBO'; %Sine with dwell, first turn to the right, 80kph

vboload
%  Channel 1  = satellites
%  Channel 2  = time
%  Channel 3  = latitude
%  Channel 4  = longitude
%  Channel 5  = velocity kmh
%  Channel 6  = heading
%  Channel 7  = height
%  Channel 8  = vertical velocity kmh
%  Channel 9  = steerang
%  Channel 10 = vxcorr
%  Channel 11 = slipcorr
%  Channel 12 = event 1 time
%  Channel 13 = rms_hpos
%  Channel 14 = rms_vpos
%  Channel 15 = rms_hvel
%  Channel 16 = rms_vvel
%  Channel 17 = latitude_raw
%  Channel 18 = longitude_raw
%  Channel 19 = speed_raw
%  Channel 20 = heading_raw
%  Channel 21 = height_raw
%  Channel 22 = vertical_velocity_raw
%  Channel 23 = true_head
%  Channel 24 = slip_angle 
%  Channel 25 = pitch_ang. 
%  Channel 26 = lat._vel.
%  Channel 27 = yaw_rate
%  Channel 28 = roll_angle 
%  Channel 29 = lng._vel.
%  Channel 30 = slip_cog
%  Channel 31 = slip_fl
%  Channel 32 = slip_fr
%  Channel 33 = slip_rl
%  Channel 34 = slip_rr
%  Channel 35 = yawrate
%  Channel 36 = x_accel
%  Channel 37 = y_accel
%  Channel 38 = temp
%  Channel 39 = pitchrate
%  Channel 40 = rollrate
%  Channel 41 = z_accel

%-----------------------------------
% SET VEHICLE DATA FOR THE VOLVO V40
%-----------------------------------
Rt=0.312;           % Tyre radius (m)
lf=0.41*2.55;       % Distance from CoG to front axis (m)
lr=2.55-lf;         % Distance from CoG to rear axis (m)
L=lf+lr;            % Wheel base (m)
h=0.2*L;            % Hight from ground to CoG (m)
mass=1435-80;       % Mass (kg)
Iz=2380;            % Yaw inertia (kg-m2)
tw=1.565;           % Track width (m)
Ratio=17;           % Steering gear ratio
Cf=80000;          % Lateral stiffness front axle (N/rad) [FREE TO TUNE]
Cr=85000;          % Lateral stiffness rear axle (N/rad) [FREE TO TUNE]
Lx_relax=0.05;      % Longitudinal relaxation lenth of tyre (m)
Ly_relax=0.15;      % Lateral relaxation lenth of tyre (m)
Roll_res=0.01;      % Rolling resistance of tyre
rollGrad=5*(pi/180);% Rollangle rad per g (rad/g)
rx=0.4;             % Distance from IMU to CoG x-axle (m)
ry=0;               % Distance from IMU to CoG y-axle (m)
rz=0;               % Distance from IMU to CoG z-axle (m)

%--------------------------------------
% SET ENVIRONEMNTAL PARAMETERS FOR TEST
%--------------------------------------
Mu=0.95;             % Coefficient of friction
g=9.81;             % Gravity constant (m/s^2)


%--------------------------------------------
% SET VARIABLES DATA FROM DATA READ FROM FILE
%--------------------------------------------
% trim_start=1;
% trim_end=length(vbo.channels(1, 2).data);
% 
% Time=(vbo.channels(1, 2).data(trim_start:trim_end,1) - vbo.channels(1, 2).data(1,1));
% yawRate_VBOX = vbo.channels(1, 35).data(trim_start:trim_end,1).*(-pi/180); %signal is inverted hence (-)
% vx_VBOX = vbo.channels(1, 5).data(trim_start:trim_end,1)./3.6;
% vy_VBOX = vbo.channels(1, 26).data(trim_start:trim_end,1)./3.6;
% ax_VBOX = vbo.channels(1, 36).data(trim_start:trim_end,1).*g;
% ay_VBOX = vbo.channels(1, 37).data(trim_start:trim_end,1).*g;
% Beta_VBOX = vbo.channels(1, 30).data(trim_start:trim_end,1).*(pi/180);
% SWA_VBOX=vbo.channels(1, 9).data(trim_start:trim_end,1).*(pi/180)/Ratio;

% % Taking away spikes in the data
% for i=1:length(Time)
%     if (i>1)
%         if (abs(SWA_VBOX(i,1)-SWA_VBOX(i-1))>1 || abs(SWA_VBOX(i,1))>7)
%             SWA_VBOX(i,1)=SWA_VBOX(i-1);
%         end
%     end
% end
n = length(Time);
dt = Time(2)-Time(1);

%----------------------------------------------
% SET MEASUREMENT AND PROCESS NOICE COVARIANCES
%----------------------------------------------
% Use as starting value 0.1 for each of the states in Q matrix
Q=(0.1)*eye(3)
% Q=[1 0 0; 0 0.025 0;0 0 0.45];%eye(3)*(0.01);%process noise
% Use as starting value 0.01 for each of the measurements in R matrix
% R=[5000 0 0;0 0.001 0; 0 0 10];
R = eye(3)*(0.01);%measurement noise 

%--------------------------------------------------
% SET INITIAL STATE AND STATE ESTIMATION COVARIANCE
%--------------------------------------------------
x_0 = [vx_VBOX(1); vy_VBOX(1); yawRate_VBOX(1)];
P_0 = Q; 


%-----------------------
% INITIALISING VARIABLES
%-----------------------
% a=eye(3);
M = x_0;
y = [vx_VBOX';ay_VBOX';yawRate_VBOX'];
% alpha=[0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000];
% beta1=[0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000];
alpha=0.04;
beta1=0.25;
% alpha=[0.0000001];
% beta1=[0.0000001];
kappa=0;%10^(-8);
%Parameters that might be needed in the measurement and state functions are added to predictParam
predictParam.dt=dt; 
updateParam.dt=dt;
% predictParam.Cf=Cf;
% predictParam.Cr=Cr;
% mat=0;

% Handles to state and measurement model functions.
state_func_UKF = @Vehicle_state_eq;
meas_func_UKF = @Vehicle_measure_eq;

a=state_func_UKF;
h=meas_func_UKF;
%-----------------------
% FILTERING LOOP FOR UKF 
%-----------------------
disp(' ');
disp('Filtering the signal with UKF...');
j=1;
k=1;
% for j=1:8
%     for k=1:8
        for i = 2:n
   
    delta=SWA_VBOX(i);
%     predictParam.delta=delta;
%     updateParam.delta=SWA_VBOX;
%     xx=cov(M(2,i),vy_VBOX(i))
%     Q=[0.1 0 0;0 0.1*cov(M(2,i),vy_VBOX(i)) 0;0 0 0.1*cov(M(3,i),yawRate(i))]%process noise
% %     Use as starting value 0.01 for each of the measurements in R matrix
%     R=(0.01).*eye(3)%measurement noise 
%     P_0=Q; 
%ukf predict
% [M(:,i),P_0] = ukf_predict1(M(:,i-1),P_0,a,Q,predictParam,alpha(j),beta1(k),kappa);
[M(:,i),P_0] = ukf_predict1(M(:,i-1),P_0,a,Q);
%ukf update

% [M(:,i),P_0,K,MU,S,LH] = ukf_update1(M(:,i),P_0,y(:,i),h,R,updateParam,alpha(j),beta1(k),kappa);
[M(:,i),P_0,K,MU,S,LH] = ukf_update1(M(:,i),P_0,y(:,i),h,R);
    % ad your predict and update functions, see the scripts ukf_predict1.m
    % and ukf_update1.m

%     
%     if i==round(n/4)
%         disp(' ');
%         disp('1/4 of the filtering done...');
%         disp(' ');
%     end
%     if i==round(n/2)
%         disp(' ');
%         disp('1/2 of the filtering done...');
%         disp(' ');
%     end
%     if i==round(n*(3/4))
%         disp(' ');
%         disp('3/4 of the filtering done... Stay tuned for the results...');
%         disp(' ');
%     end
end

%----------------------------------------
% CALCULATE THE SLIP ANGLE OF THE VEHICLE
%----------------------------------------

betaUKF(1,:)=atan(M(2,:)./M(1,:));

% for i=1:3317
%     beta1(i)=atan(M(2,i)./M(1,i));
% end
%---------------------------------------------------------
% CALCULATE THE ERROR VALES FOR THE ESTIMATE OF SLIP ANGLE
%---------------------------------------------------------
Beta_VBOX_smooth=smooth(Beta_VBOX,0.01,'rlowess'); 
[e_beta_mean,e_beta_max,time_at_max,error] = errorCalc(betaUKF',Beta_VBOX_smooth);

% disp(' ');
% fprintf('The MSE of Beta estimation is: %d \n',e_beta_mean);
% fprintf('The Max error of Beta estimation is: %d \n',e_beta_max);
% mse(j,k)=e_beta_mean;
%     end 
% end

%-----------------
% PLOT THE RESULTS
%-----------------
figure(1);
plot(Time, Beta_VBOX_smooth, Time, betaUKF);
hold on
% plot(Time,M(1,:))
% hold on
% plot(Time,vx_VBOX)
% figure (2);
% plot(Time,M(2,:))
% legend('vy_ukf')
% hold on
% plot (Time,vy_VBOX)
% legend('vy_vbox')