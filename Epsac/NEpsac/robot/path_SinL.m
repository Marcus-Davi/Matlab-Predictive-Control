% Amplitude da Senoidal, distância X da senoidal, distancia trecho horz, distancia trecho vert, v0, ts
function [Xr,Ur,iterations] = path_SinL(ampl,sin_dist,dist1,dist2,v0,Ts,x0)

X = x0;
% w0 = v0/R;


sine_iterations = round(sin_dist/(v0*Ts))
line_iterations1 = round((dist1/v0)/Ts);
line_iterations2 = round((dist2/v0)/Ts);
Tsim = (sine_iterations + line_iterations1 + line_iterations2)*Ts;
X = x0;

iterations= round((Tsim/Ts));


x = x0;
x(3) = pi/4;
x_old = 0;
y_old = 0;
dx_old = 0;
dy_old = 0;

for k = 1:iterations  % 1 a tempo de simula��o/ tempo de integra��o   



if k > sine_iterations


if(k > (sine_iterations + line_iterations1)) %y
x(1) =x(1);
x(2) = x(2) + v0*Ts; 
x(3) = pi/2;   
else %x
x(1) = x(1) + v0*Ts;
x(2) = x(2);
x(3) = 0;  
end

    v(k) = v0;
    w(k) = 0;
else %sine

x(1) = x(1) + v0*Ts;
x(2) = x(2) + (3*pi/4)*(ampl/sin_dist)*(v0*Ts)*sin(6*pi*x(1)/(4*sin_dist) + pi/2);


end

dx = x(1) - x_old;
dy = x(2) - y_old;


x_old = x(1);
y_old = x(2); 


ddx = dx - dx_old;
ddy = dy - dy_old;

dx_old = dx;
dy_old = dy;


x(3) = atan2(dy,dx);
v(k) = sqrt(dx*dx + dy*dy)/Ts;
w(k) = ((dx*ddy - dy*ddx)/(dx*dx + dy*dy))/Ts;
if(w(k) > 5)
    w(k) = 0;
end
    
X = [X x];
    

    
end
X(:,1) = X(:,2);
% x_k = X(1,:); 
% dx = x_k - [0 x_k(1:end-1)];
% 
% y_k = X(2,:); 
% dy = y_k - [0 y_k(1:end-1)];
% 
% t_k = atan2(dy,dx);
% 
% t_k(1) = t_k(2);


Xr = X(:,1:k);
% Xr(3,:) = t_k(1,1:k);
Ur = [v;w];

end