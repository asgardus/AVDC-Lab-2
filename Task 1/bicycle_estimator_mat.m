function xdot=bicycle_estimator(t,x)

global vbox_file_name Cf Cr lf lr mass Iz WA_VBOX vx_VBOX Time

r=find(Time>=t);%,1,'first');
delta_t = WA_VBOX(r(1));
vx = vx_VBOX(r(1));

xdot(1,:)=x(2); % y-speed
xdot(2,:)=(((-Cf*(((x(2)+x(4)*lf)/vx)-delta_t))-(Cr*(x(2)-x(4)*lr)/vx))/mass)-(x(4)*vx); % y-acceleration
xdot(3,:)=x(4); % Yaw rate angular velocity
xdot(4,:)=((-Cf*lf*(((x(2)+lf*x(4))/vx)-delta_t))+(Cr*lr*(x(2)-lr*x(4))/vx))/Iz; % Yaw rate angular acceleration