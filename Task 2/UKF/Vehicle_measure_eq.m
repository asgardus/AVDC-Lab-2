function y_n = Vehicle_measure_eq(x,param)
% ADDME Measurement function
%    x = the states
%    param = parameters that you might need, such as vehicle parameters.

global lf lr mass Cf Cr delta
vx=x(1,:);
vy=x(2,:);
psi_dot=x(3,:);

alpha_34=atan((vy-psi_dot*lr)/(vx));
alpha_12=atan((vy+psi_dot*lf)/(vx))-delta;

F_12=-Cf*alpha_12;
F_34=-Cr*alpha_34;

a_y=(F_34+F_12*cos(delta))/(mass);
% y_n =[vx; a_y(1,:); psi_dot];
y_n =[vx a_y psi_dot]';
%consider a_y(1,:)
