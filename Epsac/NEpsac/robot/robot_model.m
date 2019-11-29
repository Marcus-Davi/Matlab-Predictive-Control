function xk = robot_model(x0,u,Ts)
%xk -> [x,y,cos(theta)]

xk(1) = x0(1) + u(1)*cos(x0(3))*Ts;
xk(2) = x0(2) + u(1)*sin(x0(3))*Ts;
xk(3) = x0(3) + u(2)*Ts;

    % Put theta in range [-pi,pi]
    N = xk(3)/(2*pi);
 
        if N > 1;
            xk(3) = xk(3) - (round(N))*2*pi;
        end
    
        if N < -1;
            xk(3) = xk(3) - (round(N))*2*pi;
        end
    
        if xk(3) > pi
            xk(3) = xk(3) - 2*pi;
        end
    
        if xk(3) < -pi
            xk(3) = xk(3) + 2*pi;
        end
        


xk = xk';


end