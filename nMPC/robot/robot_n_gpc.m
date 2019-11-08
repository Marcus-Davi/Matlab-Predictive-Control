clear;close all;clc
%% SIMULATION PARAMETERS
Ts = 0.1;
v = 0.1;
% R = 1;
% [Xr,Ur,Tsim] = path_oito(R,v,Ts);
[Xr,Ur,Tsim] = path_S(1,v,Ts,[0 0 0]');
iterations = round((Tsim/Ts));

x0 = [0; -0.5; pi/2];
yk = x0;
vmax = 1;vmin = 0;
wmax = 1;wmin = -1;
time = 0:Ts:(Tsim-Ts);
%% NL GPC
n = 5;
nu = n;
n_in = 2;
n_out = 3;

Ql = 0.001*eye(nu*n_in); %Control

Qt = 0.01;
Q = diag([5 8 Qt]);
Qcell = repmat({Q},1,n);
Qgpc = blkdiag(Qcell{:}); %Ref
% return

%% LQR
ur1 = 0.1;
ur2 = 0.0;
A = [0 ur2 0;-ur2 0 ur1;0 0 0];
B = [1 0;0 0;0 1];
C = eye(3);
D = zeros(3,2);
SYS = ss(A,B,C,D);
Q = diag([1 100 0]);
R = eye(2);
[K_LQ,S,E] = lqr(SYS,Q,R);


%% Simula
YK = [];
UK = [];
EK = [];
uk = [0 0]';
pert = [0 0];
for k=1:(Tsim/Ts)
     
    time_s = k*Ts;
        if(time_s  > 10 )
         pert = [0 0];
        end
    
    ykest = robot_model(yk,uk,Ts); %MODELO
    yk = robot_model(yk,uk+pert,Ts); %PLANTA
    
    ne = (yk - ykest);%n(t) = y(t) - x(t)
    ni = repmat(ne,1,n); %n(t+j) = n(t+j-1)
    Yb = [];
    yb = yk;
    for i=1:n
     yb = robot_model(yb,uk,Ts) + ni(:,i);  
     Yb = [Yb; yb];
    end
    
    
        IC.x0 = yk;
        IC.u0 = uk;
      G = get_G(IC,@robot_model,1e-5,n,nu,Ts);
      
      w = [];
      uref = [];
      for m=1:n
           if(k+m > (Tsim/Ts))
              m = 0;
           end
         w = [w;Xr(:,k+m)];
         uref = [uref;Ur(:,k+m)];    
         
       
      end
%       errt = atan2(sin(errt), cos(errt)); % smooth the error

%     E = w-Yb;
%       for m=1:n
%          E( 
%       end
      
      u0 = repmat(uk,n,1);
      
      
      K0 = G'*Qgpc*G+Ql;
      K1 = G'*Qgpc*(w-Yb) + Ql*(uref-u0); %incluir velocidade
      
     duk = inv(K0)*K1;
    
     uk = uk + duk(1:n_in);
     
     ek = Xr(:,k)-yk;
    
    %thiago escolheu o erro concatenando tudo pra baixo)
    
    
    % LQR ---- INIT
    
%     e_x = Xr(1,k)-yk(1);
%     e_y = Xr(2,k)-yk(2);
%     e_t = Xr(3,k)-yk(3);
%     
%     e1 = cos(yk(3))*e_x+sin(yk(3))*e_y;
%     e2 = -sin(yk(3))*e_x+cos(yk(3))*e_y;
%     e3 = e_t;
% 
%    V = -K_LQ*[e1 e2 e3]';
%     v1 = V(1);
%     v2 = V(2);
    
%     uk(1) = Ur(1,k)*cos(e3) - v1;
%     uk(2) = Ur(2,k) -  v2;
    
    % LQR ---- END
    
    % saturation
    uk(1) = min(uk(1),vmax);
    uk(1) = max(uk(1),vmin);
    uk(2) = min(uk(2),wmax);
    uk(2) = max(uk(2),wmin);
    
    
    
    YK = [YK yk];
    UK = [UK uk];
    EK = [EK ek];
end

plot(Xr(1,:),Xr(2,:),'black--');hold on
grid on;
plot(YK(1,:),YK(2,:),'-');
title('Controle de Robô Ñ-Holonômico em trajetória')
legend('Trajetória Referência','Robô')
% plot(time,YK)
grid on;
figure
plot(time,UK);
figure
plot(time,EK);
legend('ex','ey','ez')



















function [Xr,Ur] = make_traj(Tsim,Ts)
xa = [0.5; 0.5; 0];  % x,y,theta referenciais para in�cio de trajet�ria. 
X1 = xa;         %Euler

xb = [0.5;0.5;0.0];
X2 = xb;

iterations = Tsim/Ts;

for k = 1:iterations  % 1 a tempo de simula��o/ tempo de integra��o   
    
%     %Oito
%     if(k<round(Tsim/h))
%        
%     w(k)=0.2;
%     v(k)=0.1;
%     else
%     w(k)=-0.2;
%     v(k)=0.1;
%         
%     end
    
    % "S"
    w(k)=0;
    v(k)=0.1;
    if(k<iterations/3)
        xb(3)=0;
    elseif(k<2*iterations/3 && k>=iterations/3)
        xb(3)=pi/2;
    else
        xb(3)=0;
    end



    X1 = [X1, X1(:,k)+ Ts*([v(k)*cos(xa(3));v(k)*sin(xa(3));w(k)])];  %Euler
    xa = X1(:,k+1);
    X2 = [X2, X2(:,k)+ Ts*([v(k)*cos(xb(3));v(k)*sin(xb(3));0])];  %Marcus
   % xb = X2(:,k+1);
    X2(3,k)=xb(3);
    
end


Xr = X2;
Ur = [v;w];

end







