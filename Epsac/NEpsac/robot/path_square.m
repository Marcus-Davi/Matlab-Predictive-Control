%Raio,velocidade linear (intervalo),Ts
function [Xr,Ur,Tsim] = path_square(dist,v0,Ts,x0)

% w0 = v0/R;
v_ff = linspace(v0(1),v0(2),4);

Tsim = dist/v_ff(1) + dist/v_ff(2) + dist/v_ff(3) + dist/v_ff(4);

X = x0;


iterations= round((Tsim/Ts));

itr_1 = round(dist/v_ff(1)/Ts);
itr_2 = round(dist/v_ff(2)/Ts);
itr_3 = round(dist/v_ff(3)/Ts);
% itr_4 = (dist/v_ff(4)/Ts);

for k = 1:iterations  % 1 a tempo de simula��o/ tempo de integra��o   
  
    %    "Square"
        w(k)=0;
       
        if(k<itr_1)
            xb(3)=x0(3);
             v(k)=v_ff(1);
        elseif(k>=itr_1 && k<(itr_1 + itr_2))
            xb(3)=x0(3)+pi/2;
            v(k)=v_ff(2);
        elseif(k>=(itr_1 + itr_2) && k<(itr_1 + itr_2+itr_3))
            xb(3)=x0(3)+pi;
            v(k)=v_ff(3);
        else 
            xb(3)=x0(3)+3*pi/2;
            v(k)=v_ff(4);
        end
    
    
    
    X = [X, X(:,k)+ Ts*([v(k)*cos(xb(3));v(k)*sin(xb(3));0])];  %Marcus
    % xb = X2(:,k+1);
    X(3,k)=xb(3);
    
    
end


Xr = X;
Ur = [v;w];

end