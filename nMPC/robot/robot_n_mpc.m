clear;close all;clc
%% SIMULATION PARAMETERS
Ts = 0.1;
v = 0.1;
% R = 1;
% [Xr,Ur,Tsim] = path_oito(R,v,Ts);
[Xr,Ur,Tsim] = path_S(1,v,Ts,[0 0 0]');
iterations = round((Tsim/Ts));

x0 = [0; -0.2; 0];
yk = x0;
% plot(Xr(1,:),Xr(2,:));
% return %testing

%% OPTIMIZER PARAMETERS
options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
                        'MaxIterations', 50);
%% ROBOT PARAMETERS
vmax = 1;vmin = 0;
wmax = 1;wmin = -wmax;

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
Ql = blkdiag(Qcell{:}); %Ref

Qr = 0.001;
Rl = Qr*eye(Nu*n_in);
% du = 0.00001; %increment

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
    
    [Wr,Uref] = getRef(Xr,Ur,k+1,N); 
%     [Wr,Uref] = getRef_var(Xr,Ur,k+1,N); 
    
    nf = repmat(ne,1,N);
%     nf = ne;
    u = repmat(uk,N,1);
    uk0 = uk;
    
  X = fmincon(@(u) cost(Wr,yk,Uref,u,uk0,Ts,N,Nu,[],nf,Ql,Rl),u,[],[],[],[],[],[],[],options); %SERIE-PARALELO
    

     uk = uk + X(1:2,1); %Extrai apenas o atual
     
  
    
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


function [y,naoimporta] = cost(Wr,yk,Uref,uk,uk0,Ts,N,Nu,C,nf,Ql,Rl)


n_in = length(uk)/N;
n_out= length(yk);
M_inv = eye(Nu*n_in);
n_ones = -1*ones(n_in*(Nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
M = inv(M_inv);

Yb = zeros(N*n_out,1);
yb = yk;


uref = Uref(1:n_in*Nu);
u0 = uk(1:n_in*Nu)
% Yb = [];
    for i=1:N
        u = get_block(uk,i,1,[2 1]); %pega os u's futuros
%      yb = robot_model(yb,u,Ts*i)+ nf(:,i); %Variable Horizon
     yb = robot_model(yb,u,Ts)+ nf(:,i); % Normal
     Yb = set_block(Yb,i,1,[3 1],yb);
    end
% u = uk(1:n_in*Nu);     

% L = [Uref(1:2) - uk0 ; zeros(n_in*Nu,1)];
% Ku = M_inv*EU + L(1:n_in*Nu,:);


y = (Wr-Yb)'*Ql*(Wr-Yb) + (M_inv*u0)


naoimporta =[];
end

function e = getErr(Y,W,n)
e = Y-W;
%suavização do erro em theta
for i=1:n
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










