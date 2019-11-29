clear;close all;clc
%% SIMULATION PARAMETERS
Ts = 0.1;
% R = 1;
v = 0.15;
x0 = [0; 0; 0];
[Xr,Ur,iterations] = path_L(1.5,v,Ts,x0);
real_iterations = iterations - 5^2;
% plot(Xr(1,:),Xr(2,:));
% return %testing



%% ROBOT PARAMETERS
vmax = 0.8;vmin = -vmax;
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


Qx = 10;
Qy = 10;
Qt = .00;
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
alfax = 0.99;
alfay = alfax;
alfatheta = alfax;
        % disturbance filtering x
        D = [1 -1 zeros(1,ordemFilter-1)];
        
        if (ordemFilter == 1)
            
            Cx = [1 -alfax];
            
        elseif (ordemFilter == 2)
            Cx = conv([1 -alfax],[1 -conj(alfax)]);
            
        elseif (ordemFilter == 3)
            Cx = conv([1 -alfaReal],conv([1 -alfax],[1 -conj(alfax)]));
            
        end
        
        [~,Fx] = GetFE(Cx,D,5*N);
        Fx = [Fx(1,:); Fx(3,:); Fx(6,:); Fx(10,:); Fx(15,:)]
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
        
        [~,Fy] = GetFE(Cy,D,5*N);
        Fy = [Fy(1,:); Fy(3,:); Fy(6,:); Fy(10,:); Fy(15,:)]
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
        
        [~,Ftheta] = GetFE(Ctheta,D,5*N);
        Ftheta = [Ftheta(1,:); Ftheta(3,:); Ftheta(6,:); Ftheta(10,:); Ftheta(15,:)]
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
yk = [0 -0.5 pi/2];
ykest = yk;
ykm = x0;
ykalman = x0;

% uk = [0 0]';
pert = [0 0];
% Creates noise profile
Mean = 0; % zero mean
sd_xy = 0.5; % standard deviation
sd_t = 0.0; % standard deviation
noise_xy = Mean + sd_xy.*randn(2,iterations);
noise_t = Mean + sd_t.*randn(1,iterations);
noise = [noise_xy;noise_t];


% for i=1:iterations/2
%    noise(:,i) = [0 0 0]'; 
% end

%% Simulation

for k=1:iterations

    ykest = robot_model(ykest,uk,Ts); %MODELO 
    yk = robot_model_real(yk,uk,Ts,uncertainty); % PLANTA
    
    if(k == round(iterations/2))
        yk = yk + [0.3 0 0]';
    end
       
    ykm = yk + noise(:,k);

    n = (ykm - ykest);% +  noise(:,k) ;%n(t) = y(t) - x(t)
    ykest = ykm;
        

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
    
    % Preditctions
    ub = Uref(1:n_in,1);
    Yb = zeros(n_out*N,1);    
    yb = ykest; 
 
    for j=1:N
     yb = robot_model(yb,ub,Ts*j)+ nfiltro(:,j);
     Yb = set_block(Yb,j,1,[n_out 1],yb);
    end
    

     
    % condições inicias pro impulso
        IC.x0 = ykest;
        IC.u0 = ub;
        %Toma G do modelo.
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
     
     Ku = M*EU + L(1:n_in*Nu,:);    % mgn
     
            
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
    
    
    ek = Wr(1:n_out)-yk; %plotagem
    YK = [YK yk];
    UK = [UK uk];
     EK = [EK ek];
    YKEST= [YKEST ykest];
    YKNOISE = [YKNOISE ykm];
    

end


%% PLOTS
close all
plot(Xr(1,:),Xr(2,:),'black--');hold on
grid on;
plot(YK(1,:),YK(2,:),'blue');
plot(YKEST(1,:),YKEST(2,:),'red*');
% plot(YKNOISE(1,:),YKNOISE(2,:),'magenta');
% plot(YK_Noiseless(1,:),YK_Noiseless(2,:))
title('Controle de Robô Ñ-Holonômico em trajetória')
legend('Trajetória Referência','Real Robot','YKEST')
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


