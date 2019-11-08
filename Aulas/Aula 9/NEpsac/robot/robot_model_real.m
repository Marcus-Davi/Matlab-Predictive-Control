function xk = robot_model_real(x0,u,Ts,uncertainty)
%xk -> [x,y,cos(theta)]

v_uncertain = u(1) * uncertainty.R;
w_uncertain = u(2) * uncertainty.R/uncertainty.D;

xk(1) = x0(1) + v_uncertain*cos(x0(3))*Ts;
xk(2) = x0(2) + v_uncertain*sin(x0(3))*Ts;
xk(3) = x0(3) + w_uncertain*Ts;

xk = xk';


end