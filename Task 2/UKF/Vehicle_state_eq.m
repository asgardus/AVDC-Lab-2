function x_n = Vehicle_state_eq(x,param)
% ADDME Dynamic model function
%    x = the states
%    param = parameters that you might need, such as vehicle parameters.

global lf lr mass Iz Cf Cr dt delta

vx=x(1,:);
vy=x(2,:);
psi_dot=x(3,:);

alpha_34=atan((vy-psi_dot*lr)/vx);
alpha_12=atan((vy+psi_dot*lf)/vx)-delta;

F_12=-Cf*alpha_12;
F_34=-Cr*alpha_34;

vxdot=psi_dot.*vy-(1/mass)*(F_12*sin(delta));
vydot=-psi_dot.*vx+(1/mass)*(F_34+F_12*cos(delta));
psi_ddot=(1/Iz)*(lf*F_12*cos(delta)-lr*F_34);

f_x = [vxdot; vydot; psi_ddot];


% Integrate using Runge Kutta (in the script folder) or simple euler forward

f = @(x)[f_x(1,:);f_x(2,:);f_x(3,:)];
x_n = rk4(f,dt,x(1:3,:));