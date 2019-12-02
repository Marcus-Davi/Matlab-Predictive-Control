clear;close all;clc

%% SIMULATION PARAMETERS
Ts = 0.1;
% R = 1;
v = 0.1;
x0 = [0; 0; 0];
[Xr,Ur,iterations] = path_SinL(1,1,v,Ts,x0);

% plot(Xr(1,:),Xr(2,:));
% figure
% plot(Xr(3,:));
% figure
% plot(Ur(1,:));hold on
% plot(Ur(2,:))
% return %testing

%% ROS STUFF
pub = rospublisher('/nanook_move'); %pra mover nanook
msg = rosmessage('geometry_msgs/Twist'); % msg = rosmessage(pub);
sensors = rossubscriber('/sensors'); %pegar velocidade
r = rosrate(1/Ts);
tf_tree = rostf;


%% ROBOT PARAMETERS
vmax = 0.8;vmin = -0;
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
Nu = 1;
n_in = 2;
n_out = 3;


Qx = 1;
Qy = 1;
Qt =  0.001;
Q = diag([Qx Qy Qt]);
Qcell = repmat({Q},1,N);
Qepsac = blkdiag(Qcell{:}); %Ref

Qr = 0.001;
Repsac = Qr*eye(Nu*n_in);

du = 1e-5; %increment

M_inv = eye(Nu*n_in);
n_ones = -1*ones(n_in*(Nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
M = inv(M_inv);
%% ESPAC FILTERS
ordemFilter = 2;
alfax = 0.98;
alfay = alfax;
alfatheta = 0.98;
seq = getSequence(N);
        % disturbance filtering x
        D = [1 -1 zeros(1,ordemFilter-1)];
        
        if (ordemFilter == 1)
            
            Cx = [1 -alfax];
            
        elseif (ordemFilter == 2)
            Cx = conv([1 -alfax],[1 -conj(alfax)]);
            
        elseif (ordemFilter == 3)
            Cx = conv([1 -alfaReal],conv([1 -alfax],[1 -conj(alfax)]));
            
        end
        
        [~,Fx] = GetFE(Cx,D,N*N);
        Fx = Fx([seq],:)
        npx = zeros(length(Cx)-1,1);
        [~,zix] = filter(1,Cx,0);
        
        % disturbance filtering y
        
        if (ordemFilter == 1);
            Cy = [1 -alfay];
        elseif (ordemFilter == 2);
            Cy = conv([1 -alfay],[1 -conj(alfay)]);
        elseif (ordemFilter == 3);
            Cy = conv([1 -alfaReal],conv([1 -alfay],[1 -conj(alfay)]));
        end
        
        [~,Fy] = GetFE(Cy,D,N*N);
        Fy = Fy([seq],:)
        npy = zeros(length(Cy)-1,1);
        [~,ziy] = filter(1,Cy,0);
        
        % disturbance filtering theta
        if (ordemFilter == 1);
            Ctheta = [1 -alfatheta];
        elseif (ordemFilter == 2);
            Ctheta = conv([1 -alfatheta],[1 -conj(alfatheta)]);
        elseif (ordemFilter == 3);
            Ctheta = conv([1 -alfaReal],conv([1 -alfatheta],[1 -conj(alfatheta)]));
        end
        
        [~,Ftheta] = GetFE(Ctheta,D,N*N);
        Ftheta = Ftheta([seq],:)
        nptheta = zeros(length(Ctheta)-1,1);
        [~,zitheta] = filter(1,Ctheta,0);

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
YKALM = [];
YKEST = [];
YKNOISE = [];
yk = [0 0 0];
ykest = yk;
ykm = x0;

% uk = [0 0]';
pert = [0 0];
% Creates noise profile
Mean = 0; % zero mean
sd_xy = 0.; % standard deviation
sd_t = 0.0; % standard deviation
noise_xy = Mean + sd_xy.*randn(2,iterations);
noise_t = Mean + sd_t.*randn(1,iterations);
noise = [noise_xy;noise_t];


%% Simulation

for k=1:iterations

    % ODOM
    [vread,wread] = speedGet(sensors);
%         yk = robot_model_real(yk,[vread wread],Ts,uncertainty); % PLANTA  
%SLAM
        waitForTransform(tf_tree, 'map', 'odom');
    TF = getTransform(tf_tree, 'map', 'odom');
    yk_slam = [TF.Transform.Translation.X TF.Transform.Translation.Y 0]';
    quat = [  TF.Transform.Rotation.W TF.Transform.Rotation.X TF.Transform.Rotation.Y TF.Transform.Rotation.Z];
    eul = quat2eul(quat);
    yk_slam(3) = eul(1);
    yk = yk_slam;
    
       
    ykm = yk + noise(:,k);
    
    
    
    ykest = robot_model(ykest,[vread wread],Ts); %MODELO 
    n = (ykm - ykest);% +  noise(:,k) ;%n(t) = y(t) - x(t)
    
        

    % filter x
    [nfx,zfx] = filter(1,Cx,n(1),zix);
    zix = zfx;
    npx = [nfx;npx(1:end-1)];
    Nx = Fx*npx;
    
    % filter y
    [nfy,zfy] = filter(1,Cy,n(2),ziy);
    ziy = zfy;
    npy = [nfy;npy(1:end-1)];
    Ny = Fy*npy;
    
    % filter theta
    [nftheta,zftheta] = filter(1,Ctheta,n(3),zitheta);
    zitheta = zftheta;
    nptheta = [nftheta;nptheta(1:end-1)];
    Ntheta = Ftheta*nptheta;
    
    nfiltro = [Nx'; Ny'; Ntheta'];
    
      % Pega referencias futuras
       [Wr,Uref] = getRef_var(Xr,Ur,k+1,N); 
%        [Wr,Uref] = getRef(Xr,Ur,k+1,N); 
    
    % Preditctions
    ub = Uref(1:n_in,1);
    Yb = zeros(n_out*N,1); 
    
    %% TIPOS DE RELIMENTAÇÃO
%     --- SERIE-PARALELO BEGIN ---
    yb = ykm;  %measured
    for j=1:N
     yb = robot_model(yb,ub,Ts*j)+ nfiltro(:,j);
%      yb = robot_model(yb,ub,Ts)+ nfiltro(:,j);
     Yb = set_block(Yb,j,1,[n_out 1],yb);
    end
    ykest = ykm;
% --- SERIE-PARALELO END ---
    
    
    % --- PARALELO BEGIN ---
%      yb = ykest;
%      for j=1:N 
%      yb = robot_model(yb,ub,Ts*j);
%      Yb = set_block(Yb,j,1,[n_out 1],yb);
%      end
% %     Yb = Yb  + repmat(n,N,1);
%   Yb = Yb  + reshape(nfiltro,N*n_out,1);
%     --- PARALELO END ---
    

%%      
    % condições inicias pro impulso
        IC.x0 = ykest;
        IC.u0 = ub;
        %Toma G do modelo.
%       G = get_G(IC,@robot_model,du,repmat([0 0 0]',1,N),N,Nu,Ts);
%       G = get_G_var(IC,@robot_model,du,repmat([0 0 0]',1,N),N,Nu,Ts);
%       G = get_G(IC,@robot_model,du,nfiltro,N,Nu,Ts);
      G = get_G_var(IC,@robot_model,du,nfiltro,N,Nu,Ts);

        

     
     %Calcula o erro da trajetória e base
     E = getErr(Yb,Wr,N);
     Ub = repmat(ub,Nu,1);
     
     %Erro de velocidades de referência e base
     EU = (Ub-Uref(1:2*Nu));
%       UbUr = xbnep-reshape(uref(:,1:Nu),[2*Nu,1]);
    if k>1
     L = [Ur(:,k)-uk; zeros(2*Nu,1)];
    else
    L = [-uk; zeros(2*Nu,1)];
    end
     
     Ku = M_inv*EU + L(1:n_in*Nu,:);    % mgn
     
            
      K0 = 2*(G'*Qepsac*G + M_inv'*Repsac*M_inv);
      K1 = 2*(G'*Qepsac*E + M_inv'*Repsac*Ku); %incluir velocidade
      
     Uo = -K0\K1; %Solução Analítica
    
     Ub = Ub + Uo;
     
     uk = Ub(1:n_in);
    
     
  
    
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
    
    
    
%     motorGo(pub,uk(1),uk(2))
    
    ek = Wr(1:n_out)-yk; %plotagem
    YK = [YK yk];
    UK = [UK uk];
    EK = [EK -ek];
    YKEST= [YKEST ykest];
    YKNOISE = [YKNOISE ykm];
    
    waitfor(r)
    r.statistics
    

end
  motorGo(pub,0.0,0)

%% PLOTS
close all
plot(YKNOISE(1,:),YKNOISE(2,:),'green o');
hold on
plot(Xr(1,:),Xr(2,:),'black--');
plot(YK(1,:),YK(2,:),'blue');
xlabel('X (m)','interpreter','latex')
ylabel('Y (m)','interpreter','latex')


% plot(YKNOISE(1,:),YKNOISE(2,:),'magenta');
% plot(YK_Noiseless(1,:),YK_Noiseless(2,:))
title('Controle de Robô Ñ-Holonômico em trajetória')
legend('Measured','Trajetória Referência','Real Robot')
% plot(time,YK)
time = 1:iterations;
grid on;
figure
subplot(2,1,1) %V
plot(time*Ts,UK(1,:)); hold on; 
plot(time*Ts,Ur(1,:),'black --');
ylabel('$v\ (m/s)$','interpreter','latex')
subplot(2,1,2) %omega
plot(time*Ts,UK(2,:)); hold on;
plot(time*Ts,Ur(2,:),'black --');
ylabel('$\omega\ (rad/s)$','interpreter','latex')
xlabel('Time(s)','interpreter','latex')

figure
subplot(3,1,1)
plot(time*Ts,EK(1,:));
ylabel('$(x_r - x) (m)$','interpreter','latex')
subplot(3,1,2)
plot(time*Ts,EK(2,:)); 
ylabel('$(y_r - y) (m)$','interpreter','latex')
subplot(3,1,3)
plot(time*Ts,EK(3,:)); 
ylabel('$(\theta_r - \theta) (m)$','interpreter','latex')
xlabel('Time(s)','interpreter','latex')

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


