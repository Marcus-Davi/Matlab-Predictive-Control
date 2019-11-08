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
% plot(Xr(1,:),Xr(2,:));
% return %testing

%% OPTIMIZER PARAMETERS
options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
                        'MaxIterations', 200);
%% ROBOT PARAMETERS
vmax = 0.8;vmin = 0;
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
N = 3;
Nu = 3;
n_in = 2;
n_out = 3;

Ql = 0.01*eye(Nu*n_in); %Control

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
%% ESPAC FILTERS

% % disturbance filtering x
% alfax = 0.0;
% Dx = [1 -1];
% Cx = conv([1 -alfax],[1 -alfax]);
% Dx = [Dx 0];
% [Ex,Fx] = GetFE(Cx,Dx,N);
% npx = zeros(length(Cx)-1,1);
% [xx,zix] = filter(1,Cx,0);
% 
% % disturbance filtering y
% alfay = alfax;
% Dy = [1 -1];
% Cy = conv([1 -alfay],[1 -alfay]);
% Dy = [Dy 0];
% [Ey,Fy] = GetFE(Cy,Dy,N);
% npy = zeros(length(Cy)-1,1);
% [xx,ziy] = filter(1,Cy,0);
% 
% % disturbance filtering theta
% alfatheta = 0.0;
% Dtheta = [1 -1];
% Ctheta = conv([1 -alfatheta],[1 -alfatheta]);
% Dtheta = [Dtheta 0];
% [Etheta,Ftheta] = GetFE(Ctheta,Dtheta,N);
% nptheta = zeros(length(Ctheta)-1,1);
% [xx,zitheta] = filter(1,Ctheta,0);
% 
% % return

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
% uk = [0 0]';
pert = [0 0];

% Creates noise profile
Mean = 0; % zero mean
sd_xy = 0.; % standard deviation
sd_t = 0.000; % standard deviation
noise_xy = Mean + sd_xy.*randn(2,iterations);
noise_t = Mean + sd_t.*randn(1,iterations);
noise = [noise_xy;noise_t];


for i=1:iterations/2
   noise(:,i) = [0 0 0]'; 
end

for k=1:iterations
     
    time_s = k*Ts;
    
    ykest = robot_model(yk,uk,Ts); %MODELO
    yk = robot_model_real(yk,uk,Ts,uncertainty); %PLANTA
    
    ykm = yk + noise(:,k); %yk measured
    
    ne = (ykm - ykest);% +  noise(:,k) ;%n(t) = y(t) - x(t)
    
    [Wr,Uref] = getRef(Xr,Ur,k,Nu); 
    
    nf = repmat(ne,1,N);
%     nf = ne;
    u = repmat(uk,1,Nu);
    
    X = fmincon(@(u) cost(Wr,yk,nf,u,uk,Ts,N,Nu,C,Ql),u,[],[],[],[],[],[],[],options); %SERIE-PARALELO
    

     uk = X(:,1); %Extrai apenas o atual
     
  
    
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
%     
%     uk(1) = Ur(1,k)*cos(e3) - v1;
%     uk(2) = Ur(2,k) -  v2;
    
    % LQR ---- END
    
    % saturation
    uk(1) = min(uk(1),vmax);
    uk(1) = max(uk(1),vmin);
    uk(2) = min(uk(2),wmax);
    uk(2) = max(uk(2),wmin);
    
    
    ek = Xr(:,k)-yk; %plotagem
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
time = 1:iterations;
grid on;
figure
plot(time*Ts,UK);
figure
plot(time*Ts,EK);
legend('ex','ey','ez')






%% Funções Auxiliares


function [y,naoimporta] = cost(ref,yk,nf,u,uk,Ts,n,nu,C,Ql)
x = zeros(n,1);

[~,n_in] = size(u);

M_inv = eye(nu*n_in);
n_ones = -1*ones(n_in*(nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;

nlin = length(M_inv);

u0 = [uk; zeros(nlin-1,1)];

Yb = zeros(9,1);
yb = yk;
% Yb = [];
    for i=1:n
     yb = robot_model(yb,u(:,i),Ts) + nf(i);   %incluir erro
     Yb = set_block(Yb,i,1,[3 1],yb);
% Yb = [Yb ;yb];
    end
    
y = norm(ref-Yb)^2; %+ (M_inv*u-u0)'*Ql*(M_inv*u-u0);


naoimporta =[];
end

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










