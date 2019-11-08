%Raio,velocidade linear,Ts
function [Xr,Ur,Tsim] = path_oito(R,v0,Ts,x0)

X = x0;

d = 2*pi*R*2; %duas circunferencias
w0 = v0/R;

Tsim = d/(v0);


iterations= round((Tsim/Ts));

for k = 1:iterations  % 1 a tempo de simula��o/ tempo de integra��o   
    
%     %Oito
    if(k<iterations/2)      
    w(k)=w0;
    v(k)=v0;
    else
    w(k)=-w0;
    v(k)=v0;        
    end   
    % "S"
%     w(k)=0;
%     v(k)=0.1;
%     if(k<iterations/3)
%         forced_angle=0;
%     elseif(k<2*iterations/3 && k>=iterations/3)
%         forced_angle=pi/2;
%     else
%         forced_angle=0;
%     end

    
    X = [X, X(:,k)+ Ts*([v(k)*cos(X(3,k));v(k)*sin(X(3,k));w(k)])];  
 
    
%     X = [X, X(:,k)+ Ts*([v(k)*cos(forced_angle);v(k)*sin(forced_angle);0])]; 
%     X(3,end) = forced_angle;
    
    
end


Xr = X;
Ur = [v;w];

end