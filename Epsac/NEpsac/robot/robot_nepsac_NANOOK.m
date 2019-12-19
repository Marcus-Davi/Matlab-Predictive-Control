clear;close all;clc
%% SIMULATION / EXPERIMENT PARAMETERS 
HIL = 1; %Hard-In-Loop

Ts = 0.1;
% R = 1;
v = 0.1;
x0 = [0.; 0; 0];
%4 3.5 1.2 2 -> grande
% [Xr,Ur,iterations] = path_SinL(4,3.5,1.2,2,v,Ts,x0);   % grande
[Xr,Ur,iterations] = path_SinL(2,2.5,1.2,2,v,Ts,x0);   % pequena
% [Xr,Ur,iterations] = path_square(2,[v v],Ts,x0);   
% [Xr,Ur,iterations] = path_L(2.25,v,Ts,x0);

plot(Xr(1,:),Xr(2,:));

grid on
% figure
% plot(Xr(3,:));
% figure
% plot(Ur(1,:));hold on
% plot(Ur(2,:))
% return %testing

%% ROS STUFF
if(HIL)
pub = rospublisher('/nanook_move'); %pra mover nanook
msg = rosmessage('geometry_msgs/Twist'); % msg = rosmessage(pub);
sensors = rossubscriber('/sensors'); %pegar velocidade
slam = rossubscriber('/slam_out_pose');
r = rosrate(1/Ts);
% load('MagCalibration.mat'); %MagOff
end

%% ROBOT PARAMETERS
vmax = 0.25;vmin = 0;
wmax = 0.5;wmin = -wmax;

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
du = 1e-6; %increment
M_inv = eye(Nu*n_in);
n_ones = -1*ones(n_in*(Nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
M = inv(M_inv);
%% ESPAC FILTERS
ordemFilter = 2;
alfax = 0.9;
alfay = alfax;
alfatheta = 0.9;
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

%% Simula
uk = [0 0]';
YKALM = [];
YKEST = [];
YKNOISE = [];
yk = [0 0 0]'; %real
ykest = yk; %model
yk_odom = yk;
ykm = x0; 
yb = yk;

% uk = [0 0]';
pert = [0 0];
% Creates noise profile
Mean = 0; % zero mean
sd_xy = 0.0; % standard deviation
sd_t = 0.0; % standard deviation
noise_xy = Mean + sd_xy.*randn(2,iterations);
noise_t = Mean + sd_t.*randn(1,iterations);
noise = [noise_xy;noise_t];

%% Graphics
YKM = zeros(n_out,iterations);
UK = zeros(n_in,iterations);
EK = zeros(n_out,iterations);
% return
%% Simulation
yaw_zero_angle = 0;

for k=1:iterations
    
    
    if(HIL) %Experiment
        
    % SENSOR READING
    sensor_data = receive(sensors);
    [vread,wread] = speedGet(sensor_data);
%     yaw = angleMagGet(sensor_data,MagOff);
%     if(yaw_zero_angle == 0) %magnetometer start
%         yaw_zero_angle = yaw;
%     end
    % Odometry
      yk_odom = robot_model(yk_odom,[vread wread],Ts); % PLANTA 
%       yk_odom(3) = -(yaw - yaw_zero_angle); %magnetometer
      
    % SLAM
    slam_pose = receive(slam); %SLAM   
    quat = [  slam_pose.Pose.Orientation.W slam_pose.Pose.Orientation.X, ...
              slam_pose.Pose.Orientation.Y, slam_pose.Pose.Orientation.Z];
    eul = quat2eul(quat);
    yk_slam = [slam_pose.Pose.Position.X slam_pose.Pose.Position.Y eul(1)]';
        
    
    % Pose Feedback. Use either slam or odometry
    ykm = yk_slam;
%     ykm = yk_odom;
    
    ykest = robot_model(yb,[vread wread],Ts); %MODEL ESTIMATION
    ykest(3) = wrapToPi(ykest(3));
    
    yk = ykm; %For Plots
    
    
    else %Simulation
    yk = robot_model_real(yk,uk,Ts,uncertainty); % PLANT
    yk(3) = wrapToPi(yk(3));
    ykest = robot_model(yb,uk,Ts); %MODEL ESTIMATION
    ykest(3) = wrapToPi(ykest(3));
    ykm = yk + noise(:,k);
    vread = uk(1);
    wread = uk(2);
    end
        
    %Error of estimation   
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
    
      % Get Future References
       [Wr,Uref] = getRef_var(Xr,Ur,k,N); %Variable Horizon
%        [Wr,Uref] = getRef(Xr,Ur,k,N); %Normal Horizon
    
    % Preditctions
    ub = Uref(1:n_in,1); %U Base
%     ub(2) = 0;
    Yb = zeros(n_out*N,1); %Y Base
    
    %% TIPOS DE RELIMENTAÇÃO
%     --- SERIE-PARALELO BEGIN ---
    yb = ykm;  %measured
    for j=1:N
     yb = robot_model(yb,ub,Ts*j)+ nfiltro(:,j); %Variable Horizon
%      yb = robot_model(yb,ub,Ts)+ nfiltro(:,j); % Normal Horizon
     Yb = set_block(Yb,j,1,[n_out 1],yb);
    end
    yb = ykm;
% --- SERIE-PARALELO END ---
    
    
    % --- PARALELO BEGIN ---
%      yb = ykest; %estimated
%      for j=1:N 
% %      yb = robot_model(yb,ub,Ts*j);
%      yb = robot_model(yb,ub,Ts);
%      Yb = set_block(Yb,j,1,[n_out 1],yb);
%      end
% %     Yb = Yb  + repmat(n,N,1);
%   Yb = Yb  + reshape(nfiltro,N*n_out,1);
%   yb = ykest;
%     --- PARALELO END ---
    

%%      
    % condições inicias pro impulso
        IC.x0 = yb;
        IC.u0 = ub;
        %Toma G do modelo.
%       G = get_G(IC,@robot_model,du,repmat([0 0 0]',1,N),N,Nu,Ts);
%       G = get_G_var(IC,@robot_model,du,repmat([0 0 0]',1,N),N,Nu,Ts);
%       G = get_G(IC,@robot_model,du,nfiltro,N,Nu,Ts);
      G = get_G_var(IC,@robot_model,du,nfiltro,N,Nu,Ts);

     if(k == 585)
    k
    end

if k == 185
k
end
     
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
      
    
    ek = Xr(:,k)-yk; %plotagem
    YKM(:,k) = yk;
%     UK(:,k) = uk;
    UK(:,k) = [vread wread]';
    ek(3) = atan2(sin(ek(3)),cos(ek(3)));
    EK(:,k) = -ek;
%     YKEST= [YKEST ykest];
%     YKNOISE = [YKNOISE ykm];

    %% Runtime
    if(HIL)
    motorGo(pub,uk(1),uk(2));
%     plot(Xr(1,:),Xr(2,:));
%     hold on;
%     plot(YKM(1,:),YKM(2,:));
%     hold off;
    
    waitfor(r)
    r.statistics
    end
    

end
if(HIL)
  motorGo(pub,0.0,0);
end

%% PLOTS
close all
% plot(YKNOISE(1,:),YKNOISE(2,:),'green o');
hold on
plot(Xr(1,:),Xr(2,:),'black--');
plot(YKM(1,:),YKM(2,:),'blue');
xlabel('X (m)','interpreter','latex')
ylabel('Y (m)','interpreter','latex')


% plot(YKNOISE(1,:),YKNOISE(2,:),'magenta');
% plot(YK_Noiseless(1,:),YK_Noiseless(2,:))
title('Controle de Robô Ñ-Holonômico em trajetória')
legend('Reference Trajectory','Robot')
% plot(time,YK)
time = 1:iterations;
grid on;
figure
subplot(2,1,1) %V
plot(time*Ts,UK(1,:)); hold on; 
% plot(time*Ts,Ur(1,:),'black --');
ylabel('$v\ (m/s)$','interpreter','latex')
subplot(2,1,2) %omega
plot(time*Ts,UK(2,:)); hold on;
% plot(time*Ts,Ur(2,:),'black --');
ylabel('$\omega\ (rad/s)$','interpreter','latex')
xlabel('Time(s)','interpreter','latex')

figure
subplot(3,1,1)
plot(time*Ts,EK(1,:));
ylabel('$(x_r - x) (m)$','interpreter','latex')
grid on;
subplot(3,1,2)
plot(time*Ts,EK(2,:)); 
grid on;
ylabel('$(y_r - y) (m)$','interpreter','latex')
subplot(3,1,3)
plot(time*Ts,EK(3,:)); 
grid on;
ylabel('$(\theta_r - \theta) (m)$','interpreter','latex')
xlabel('Time(s)','interpreter','latex')

%% Funções Auxiliares


function e = getErr(W,Y,nu)
e = W-Y;
%suavização do erro em theta
for i=1:nu
   e(i*3) = atan2(sin(e(i*3)),cos(e(i*3)));  %3 pq são 3 saídas! 
% 	if(abs(e(i*3)) > pi)
% 		if(e(i*3) > 0.0)
% 		e(i*3) = e(i*3) - 2*pi;
% 		else
% 		e(i*3) = e(i*3) + 2*pi;
%         end
%     end
	
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
   seq(i) = i*(i+1)/2;%1 3 6 10 15
end

end


