%% SIMULATION PARAMETERS
clear;close all;clc
Ts = 0.1;
% R = 10;
v = 0.1;
x0 = [0; 0; 0]; %-1.4516
% entrada est | saida est ~ sensor
% -0.3 | 2.8 ~ bussola
% -1.74 | 0.54 ~ GPS

% [Xr,Ur,Tsim] = path_oito(2,v,Ts,x0); 
% [Xr,Ur,Tsim] = path_reta(15,v,Ts,x0); 
[Xr,Ur,Tsim] = path_L(1,v,Ts,x0);

%ir de 0,0 -> (-6.5,-6.5)
iterations = round((Tsim/Ts));

% x0 = [0 ; 0 ; 0];

plot(Xr(1,:),Xr(2,:));
hold on;
% return %testing


%% EPSAC PARAMETERS
vmax = 0.4;vmin = 0;
wmax = 0.7;wmin = -wmax;

N = 5;
Nu = 1;
n_in = 2;
n_out = 3;

Repsac = 0.01*eye(Nu*n_in); %Control

Qt = 0.001;
Q = diag([1 1 Qt]);
Qcell = repmat({Q},1,N);
Qepsac = blkdiag(Qcell{:}); %Ref
du = 0.00001; %increment

M_inv = eye(Nu*n_in);
n_ones = -1*ones(n_in*(Nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
M = inv(M_inv);

%% LQR
% ur1 = 0.1;
% ur2 = 0.0;
% A = [0 ur2 0;-ur2 0 ur1;0 0 0];
% B = [1 0;0 0;0 1];
% C = eye(3);
% D = zeros(3,2);
% SYS = ss(A,B,C,D);
% Q = diag([1 100 0]);
% R = eye(2);
% [K_LQ,S,E] = lqr(SYS,Q,R);

%% Simulation Parameters
YK = [];
UK = [];
EK = [];
uk = Ur(:,1);
YKALM = [];
GPS = [];
YKEST = [];
YKM = [];
% yk = x0;
yk = [0 -0.5 pi/2]';
ykalman = yk;


% NOISE PROFILE
Ts_noise_xy = Ts;
Noise_time = iterations/2;
Decim = Ts_noise_xy/Ts;


Mean = 0; % zero mean
sd_xy = 0; % standard deviation
sd_t = 0.1; % standard deviation
noise_xy = [0;0];
noise_t = 0;

%% KALMAN FILTER PARAMETERS
Q_kal = Ts^2*eye(3); %0.0005 GPS, 0.05 ODO
Q_kal(3,3) = 0.001;
Q_kal = 0.001*Q_kal;
R_kal = 1*[sd_xy^2 0 0;0 sd_xy^2 0;0 0 sd_t^2];
Pk = zeros(3);

%% Simulation

% gps_i = 1;

for k=1:iterations
                
      yk = robot_model(yk,uk,Ts); %real robot
      
        if(k==round(iterations/2))
        yk = yk + [-0 0 0]';
        end
        
      [ykalman,Pk] = kalman_ext_predict(Ts,@robot_model,Q_kal,ykalman,uk,Pk);
      
      
      if(rem(k,Decim)==0 && k>Noise_time) % GPS DECIMATION   
        noise_xy = Mean + sd_xy.*randn(2,1)/2;
        GPS = [GPS; ykm(1) ykm(2)];
      ykm = yk + [noise_xy;0];
      [ykalman,~,Pk,~] = kalman_ext_update(Ts,@measurement_model,R_kal,ykalman,uk,Pk,ykm); 
      else
      ykm = yk + [noise_xy;noise_t];

      end
               
    ne_kalman = (ykm - ykalman);% +  noise(:,k) ;%n(t) = y(t) - x(t)
    
    nfiltro = repmat(ne_kalman,1,N); %@FILTER OVERRIDE
    
    %preditctions
    ub = [0.3 0]'; 
    Yb = [];

    yb = ykalman;
    for j=1:N
     yb = robot_model(yb,ub,j*Ts);
     Yb = [Yb; yb];
    end
       
%     else
%    yb = ykest; %ykest
%     for j=1:N
%         yb = robot_model(yb,ub,j*Ts) ;
%         Yb = [Yb; yb];     
%     end
%         Yb = Yb + repmat(ne_raw,N,1);%repmat(ne,N,1);
%     end
    
    %condições inicias
        IC.x0 = ykalman;
        IC.u0 = ub;
        %Toma G do modelo.
        
        
      G = get_G(IC,@robot_model,du,nfiltro,N,Nu,Ts);
%     G = get_G_var(IC,@robot_model,du,N,Nu,Ts);
      
      %Pega referencias futuras
     
      [Wr,Uref] = getRef_var(Xr,Ur,k,N); 
%         [Wr,Uref] = getRef(Xr,Ur,k,N); 
     
     %Calcula o erro da trajetória e base
     
     E = getErr(Wr,Yb,N);
     Ub = repmat(ub,Nu,1);
     
     %Erro de velocidades de referência e base
     EU = (Uref(1:2*Nu)-Ub);
            
     
     % REVISAR AQUI !!
%    K0 = 2*(G'*Qepsac*G+M_inv'*Repsac*M_inv);
%    K1 = 2*(G'*Qepsac*E + Repsac*EU); %incluir velocidade
   
      K0 = 2*(G'*Qepsac*G+Repsac);
      K1 = 2*(G'*Qepsac*E + Repsac*EU); %incluir velocidade
      
     Uo = K0\K1; %Solução Analítica
    
     uk = Ub + Uo;
     
     uk = uk(1:2); %Extrai apenas o atual
     
  
    
    % LQR OVERRIDE ---- INIT
%     
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
%     
%     uk(1) = Ur(1,k)*cos(e3) - v1;
%     uk(2) = Ur(2,k) -  v2;
    
    % LQR OVERRIDE ---- END
    
    % saturation
    uk(1) = min(uk(1),vmax);
    uk(1) = max(uk(1),vmin);
    uk(2) = min(uk(2),wmax);
    uk(2) = max(uk(2),wmin);
      
    ek = Xr(:,k)-yk; %plotagem
    YK = [YK yk];
    UK = [UK uk];
    EK = [EK ek];
    YKALM = [YKALM ykalman];
    YKM = [YKM ykm];
    

    
    

end




%% PLOTS
close all
figure;hold on;
% plot(YKM(1,:),YKM(2,:),'blue*');hold on
plot(GPS(:,1),GPS(:,2),'g*');
plot(Xr(1,:),Xr(2,:),'black--');
grid on;
plot(YKALM(1,:),YKALM(2,:),'blue');
plot(YK(1,:),YK(2,:),'red');

% plot(YK_Noiseless(1,:),YK_Noiseless(2,:))
title('Controle de Robô Ñ-Holonômico em trajetória')
legend('GPS','Trajetória Referência','Kalman','Real Robot')
% plot(time,YK)
time = 1:iterations;
grid on;

figure
plot(time*Ts,UK);
grid on;

figure
plot(time*Ts,EK);
legend('ex','ey','e_\theta')
grid on;

% figure
% plot(Xr(1,:))
% hold on;
% plot(YKALM(1,:))
% plot(YK(1,:))
% legend('Ref','Kalman','Real')

    




%% Funções Auxiliares

function e = getErr(W,Y,nu)
e = W-Y;
%suavização do erro em theta
for i=1:nu
   e(i*3) = atan2(sin(e(i*3)),cos(e(i*3)));  %3 pq são 3 saídas!     
end
end

function [wr,ur] = getRef(Xr,Ur,k,nu)
wr = [];
ur = [];
[~,kend] = size(Xr);
for i=1:nu    
     %testa final da traj. se sim, repete
    if(k+i < kend)
    wr = [wr;Xr(:,k+i)];    
    ur = [ur;Ur(:,k+i)];    
    else
    wr = [wr;Xr(:,end)];
    ur = [ur;Ur(:,end)];
    end
end

end

function [wr,ur] = getRef_var(Xr,Ur,k,nu)
wr = [];
ur = [];
[~,kend] = size(Xr);
Seq = getSequence(nu);
for i=1:nu    
     %testa final da traj. se sim, repete
    if(k+Seq(i) < kend)
    wr = [wr;Xr(:,k+Seq(i))];    
    ur = [ur;Ur(:,k+Seq(i))];    
    else
    wr = [wr;Xr(:,end)];
    ur = [ur;Ur(:,end)];
    end
end

end

function seq = getSequence(n)
seq = zeros(n,1);
for i=1:n
   seq(i) = i*(i+1)/2;
end

end









