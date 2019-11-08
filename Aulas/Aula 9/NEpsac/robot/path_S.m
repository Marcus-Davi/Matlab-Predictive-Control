%Raio,velocidade linear,Ts
function [Xr,Ur,Tsim] = path_S(dist,v0,Ts,x0)

X = x0;
% w0 = v0/R;
w0 = 0;

Tsim = 3*dist/(v0);

X = x0;

iterations= round((Tsim/Ts));

for k = 1:iterations  % 1 a tempo de simula��o/ tempo de integra��o   
  
    %    "Z"
        w(k)=0;
        v(k)=v0;
        if(k<iterations/3)
            xb(3)=x0(3);
        elseif(k>=iterations/3 && k<2*iterations/3)
            xb(3)=x0(3)+pi/2;
        else
            xb(3)=x0(3);
        end
    
    
    
    X = [X, X(:,k)+ Ts*([v(k)*cos(xb(3));v(k)*sin(xb(3));0])];  %Marcus
    % xb = X2(:,k+1);
    X(3,k)=xb(3);
    
    
end


Xr = X;
Ur = [v;w];

end