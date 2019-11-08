%Raio,velocidade linear,Ts
function [Xr,Ur,Tsim] = path_reta(dist,v0,Ts,x0)

X = x0;
% w0 = v0/R;
w0 = 0;

Tsim = dist/(v0);



iterations= round((Tsim/Ts));

for k = 1:iterations  % 1 a tempo de simula��o/ tempo de integra��o   
    v(k) = v0;
    w(k) = 0;


    
    X = [X, X(:,k)+ Ts*([v(k)*cos(X(3,k));v(k)*sin(X(3,k));w(k)])];  
 
    
%     X = [X, X(:,k)+ Ts*([v(k)*cos(forced_angle);v(k)*sin(forced_angle);0])]; 
%     X(3,end) = forced_angle;
    
    
end


Xr = X;
Ur = [v;w];

end