function [y] = modelo(u,xo,h)
%This function represents the robot model
%   Detailed explanation goes here

 % sampling time
Ts = 0.1;
k = 1;
%h = Ts/k; %integration time

for l = 1:k;
    
    % Euler
%     M = [cos(xo(3)) 0;sin(xo(3)) 0;0 1];
%     x = xo + h*M*u;
    
%     Runge-kutta
    k1 = h*[u(1)*cos(xo(3)); u(1)*sin(xo(3));u(2)];
    k2 = h*[u(1)*cos(xo(3)+k1(3)/2); u(1)*sin(xo(3)+k1(3)/2);u(2)];
    k3 = h*[u(1)*cos(xo(3)+k2(3)/2); u(1)*sin(xo(3)+k2(3)/2);u(2)];
    k4 = h*[u(1)*cos(xo(3)+k3(3)); u(1)*sin(xo(3)+k3(3));u(2)];
    x = xo + (1/6)*(k1+2*k2+2*k3+k4);
    
    % Put theta in range [-pi,pi]
    N = x(3)/(2*pi);
 
        if N > 1;
            x(3) = x(3) - (round(N))*2*pi;
        end
    
        if N < -1;
            x(3) = x(3) - (round(N))*2*pi;
        end
    
        if x(3) > pi
            x(3) = x(3) - 2*pi;
        end
    
        if x(3) < -pi
            x(3) = x(3) + 2*pi;
        end
        
    xo = x; % redefine xo
        
end

y = xo;

end

