clear;close all;clc
%% SIMULATION PARAMETERS
Ts = 0.1;
R = 10;
v = 0.1;
x0 = [1; 1; 0];
% [Xr,Ur,Tsim] = path_oito(R,v,Ts);
% [Xr,Ur,Tsim] = path_reta(15,v,Ts);
[Xr,Ur,Tsim] = path_S(1.2,v,Ts,x0);
iterations = round((Tsim/Ts));


% x0 = [0 ; 0 ; 0];

% plot(Xr(1,:),Xr(2,:));
% return %testing

%0 -> PARAELO | 1-> SERIE_PARALELO , 0 -> PARALELO  <----- ESCOLHA AQUI ###
SERIE_PARALELO = 1; 

%% ROBOT PARAMETERS
vmax = 0.4;vmin = 0;
wmax = 0.8;wmin = -wmax;

%Model Robot
ROBOT.R = 0.08; %wheel radius
ROBOT.D = 0.4; %Robot diameter

%Real Robot
ROBOT_REAL.R = 0.08;
ROBOT_REAL.D = 0.4;

uncertainty.R = ROBOT.R/ROBOT_REAL.R;
uncertainty.D = ROBOT.D/ROBOT_REAL.D;



%% EPSAC PARAMETERS
N = 5;
Nu = 5;
n_in = 2;
n_out = 3;

Ql = 0.001*eye(Nu*n_in); %Control

Qt = 0.01;
Q = diag([1 1 Qt]);
Qcell = repmat({Q},1,N);
Qepsac = blkdiag(Qcell{:}); %Ref
du = 0.00001; %increment

M_inv = eye(Nu*n_in);
n_ones = -1*ones(n_in*(Nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
M = inv(M_inv);
%% ESPAC FILTERS THIAGO

% % disturbance filtering x
% alfax = 0.98;
% Dx = [1 -1];
% Cx = conv([1 -alfax],1);
% Dx = [Dx ];%[Dx 0]
% [Ex,Fx] = GetFE(Cx,Dx,N);
% npx = zeros(length(Cx)-1,1);
% [~,zix] = filter(1,Cx,0);
% 
% % disturbance filtering y
% alfay = alfax;
% Dy = [1 -1];
% Cy = conv([1 -alfay],1);
% Dy = [Dy ];%[Dy 0]
% [Ey,Fy] = GetFE(Cy,Dy,N);
% npy = zeros(length(Cy)-1,1);
% [~,ziy] = filter(1,Cy,0);
% 
% % disturbance filtering theta
% alfatheta = 0.0;
% Dtheta = [1 -1];
% Ctheta = conv([1 -alfatheta],[1 -alfatheta]);
% Dtheta = [Dtheta 0];
% [Etheta,Ftheta] = GetFE(Ctheta,Dtheta,N);
% nptheta = zeros(length(Ctheta)-1,1);
% [xx,zitheta] = filter(1,Ctheta,0);

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
uk = Ur(:,1);
YKALM = [];
YKEST = [];
YKNOISE = [];
yk = [0 0 0]';
ykest = yk;
ykalman = yk;

% uk = [0 0]';
pert = [0 0];
% Creates noise profile
Mean = 0; % zero mean
sd_xy = 3; % standard deviation
sd_t = 0.01; % standard deviation
noise_xy = Mean + sd_xy.*randn(2,iterations)/2;
noise_t = Mean + sd_t.*randn(1,iterations)/2;
noise = [noise_xy;noise_t];


% for i=1:iterations/2
%    noise(:,i) = [0 0 0]'; 
% end


%% KALMAN FILTER PARAMETERS
Q_kal = 0.001*eye(3); %0.0002
Q_kal(3,3)=0.00001;
% Q_kal(3,3) = 0.0001; %bussola
R_kal = [sd_xy^2 0 0;0 sd_xy^2 0;0 0 sd_t^2];
Pk = zeros(3);

%% Simulation

for k=1:iterations
     
    time_s = k*Ts;
    
    
    if(SERIE_PARALELO)
    ykest = robot_model(yk,uk,Ts); %MODELO
    else
    ykest = robot_model(ykest,uk,Ts); %MODELO    
    end
       
    yk = robot_model_real(yk,uk,Ts,uncertainty); %PLANTA
%     YK_Noiseless = [YK_Noiseless yk];
    
    if(time_s == 260)
       yk = yk + [0 10 0]';  %chute
    end
%     ykm = yk + noise(:,k); %yk measured
    ykm = yk;
    
    
    [ykalman,~,Pk,err] = kalman_ext(Ts,@robot_model,@measurement_model,Q_kal,R_kal,ykalman,uk,Pk,ykm);
    
    ne_kalman = (ykm - ykalman);% +  noise(:,k) ;%n(t) = y(t) - x(t)
    ne_raw = ykm-ykest;
    
%     ne = err;

%         % filter x
%     [nfx(k),zfx] = filter(1,Cx,ne(1),zix);
%     zix = zfx;
%     npx = [nfx(k);npx(1:end-1)];
%     Nx = Fx*npx;
% 
%     % filter y
%     [nfy(k),zfy] = filter(1,Cy,ne(2),ziy);
%     ziy = zfy;
%     npy = [nfy(k);npy(1:end-1)];
%     Ny = Fy*npy;
%     
%         % filter theta
%     [nftheta(k),zftheta] = filter(1,Ctheta,ne(1),zitheta);
%     zitheta = zftheta;
%     nptheta = [nftheta(k);nptheta(1:end-1)];
%     Ntheta = Ftheta*nptheta;
    
%      nfiltro = [Nx'; Ny'; Ntheta']; %a linha i é a predição k+i do erro

     nfiltro = repmat(ne_kalman,1,N); %@FILTER OVERRIDE
    

    
    %preditctions
    ub = [0.3 0]'; 
%     ub = Ur(:,k);
    Yb = [];
    
    if(SERIE_PARALELO)
    yb = ykalman;
    for j=1:N
     yb = robot_model(yb,ub,Ts); %+ nfiltro(:,j);
     Yb = [Yb; yb];
    end
       
    else
   yb = ykest; %ykes
    for j=1:N
        yb = robot_model(yb,ub,Ts) ;
        Yb = [Yb; yb];     
    end
        Yb = Yb + repmat(ne_raw,N,1);%repmat(ne,N,1);
    end
    
    %condições inicias
        IC.x0 = ykalman;
        IC.u0 = ub;
        %Toma G do modelo.
      G = get_G(IC,@robot_model,du,N,Nu,Ts);
      
      %Pega referencias futuras
     
      [Wr,Uref] = getRef(Xr,Ur,k,N); 
     
     %Calcula o erro da trajetória e base
     E = getErr(Wr,Yb,N);
     Ub = repmat(ub,Nu,1);
     %Erro de velocidades de referência e base
     EU = (Uref(1:2*Nu)-Ub);
            
%    K0 = 2*(G'*Qepsac*G+M_inv'*Ql*M_inv);
%    K1 = 2*(G'*Qepsac*E + Ql*EU); %incluir velocidade
       
      K0 = 2*(G'*Qepsac*G+Ql);
      K1 = 2*(G'*Qepsac*E + Ql*EU); %incluir velocidade
      
     Uo = K0\K1; %Solução Analítica
    
     uk = Ub + Uo;
     
     uk = uk(1:2); %Extrai apenas o atual
     
  
    
    % LQR OVERRIDE ---- INIT
    
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
    YKEST= [YKEST ykest];
    YKNOISE = [YKNOISE ykm];
end


%% PLOTS
close all
plot(YKNOISE(1,:),YKNOISE(2,:),'green*');hold on
plot(Xr(1,:),Xr(2,:),'black--'); hold on
grid on;
plot(YKALM(1,:),YKALM(2,:),'blue');
plot(YK(1,:),YK(2,:),'red');

% plot(YK_Noiseless(1,:),YK_Noiseless(2,:))
title('Controle de Robô Ñ-Holonômico em trajetória')
legend('Noise','Trajetória Referência','Kalman','Real Robot')
% plot(time,YK)
time = 1:iterations;
grid on;

figure
plot(time*Ts,UK);
grid on;

figure
plot(time*Ts,EK);
legend('ex','ey','ez')
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










