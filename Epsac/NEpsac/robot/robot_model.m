function xk = robot_model(x0,u,Ts)
%xk -> [x,y,cos(theta)]

xk(1) = x0(1) + u(1)*cos(x0(3))*Ts;
xk(2) = x0(2) + u(1)*sin(x0(3))*Ts;
xk(3) = x0(3) + u(2)*Ts;

xk = xk';


end