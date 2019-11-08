clear;close all;clc
addpath('../../GPC_LIB')
Ts = 60; %60s
Tsim = Ts*30; %30 x 60s
time = 0:Ts:Tsim;


%Iterações

%% Parametros Controle
n = 3;
nu = 3;
Ql = 0.1*eye(nu);

n_in = 1;

% M_inv = eye(nu*n_in);
% n_ones = -1*ones(n_in*(nu-1),1);
% n_diag = diag(n_ones,-n_in);
% M_inv = M_inv+n_diag;
%% Linearização do Modelo
g = 9.81;
Area = 10;
a = 0.01*pi;
k1 = sqrt(2*g)*a;

h0 = 5; %Ponto de operação

A = -k1/(2*Area*sqrt(h0));
B = 1/Area;
H = 1;
D = 0;
sys_linear = ss(A,B,H,D);

Q_LQ = 0.1;
R_LQ = 1;
K_LQ = lqr(sys_linear,Q_LQ,R_LQ);




% G = [];F = [];E=[];
% siz = length(H*B);
% HB_old = [];
% for i = 1:n
%     HAB = H*A^(i-1)*B;
%     M = [HAB HB_old repmat(zeros(siz,siz),[1 n-i])];
%     G = [G;M];
%     F = [F;H*A^(i)];
%     E = [E;H*A^(i-1)*D];   
%     HB_old = [HAB HB_old];
% end

% G = G(:,1:nu*1);
% K = inv(G'*G + Ql)*G';
% %dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% % como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)
% Kr = K(1,:) * ones(n,1);
% Kf = K(1,:) * F;
% Ke = K(1,:) * E;

% return

yk = h0; %h0
uk = k1*sqrt(yk);
% uk = 0

%% Modelo
YK = [];
UK = [];

pert = 0;
ref = h0;

for k=time 
    time_s = k/Ts;
        if(time_s  > 15 )
        pert = 1;
        ref = 5;
%         elseif(time_s > 3)
%             pert = 0;
        end
        
            ykest = tank(yk,uk,Ts); % RODA A PLANTA (ESTIMAÇÃO)  
            
            yk = tank_incerto(yk,uk+pert,Ts); % RODA A PLANTA
       
        
     
        
%dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)

    ne = (yk - ykest);%n(t) = y(t) - x(t)
    ni = ne*ones(n,1); %n(t+j) = n(t+j-1)
    yb = yk;
    for i=1:n
     yb = tank(yb,uk,Ts) + ne;  %+ni(i)
     Yb(i) = yb;
    end
    
      IC.u0 = uk;
      IC.x0 = yk;
        G = get_G(IC,@tank,0.01,n,nu,Ts);
    
%     K0 = inv(G'*G+M_inv'*Ql*M_inv);
      w = ref*ones(n,1);    
     duk = inv(G'*G+Ql)*(G'*(w-Yb));
  
     duk = duk(1);
   
   
%    duk = -K_LQ*(yk-h0); %estado linearizado LQR
   uk = uk+duk;

   YK = [YK yk];
   UK = [UK uk];
end
time = time/Ts;
plot(time,YK)
hold on
plot(time,UK)
grid on;
xlabel('Tempo [m]')
ylabel('Altura Tanque [m]')


%% Cost

function [y,naoimporta] = cost(ref,yk,u,Ts,n)
x = zeros(n,1);

yb = yk;
    for i=1:n
     yb = tank(yb,u(i),Ts);   %incluir erro
     x(i) = yb;
    end
    
    y = norm(ref-x)^2 ;%+ norm(u)^2;
    
naoimporta =[];
end


