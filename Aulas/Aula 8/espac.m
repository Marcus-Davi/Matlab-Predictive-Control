clear;close all;clc
%% Planta
Ts = 1;
delta = [1 -1];
Den = [1 -0.8 0];
Num = [0 0.4 0.6] ;
Pz = tf(Num,Den,Ts);
%Horizontes de Controle
n = 3;
nu = 3;
GPC = gpc_tf2ss_U(Pz,[],delta,n,nu); %Apenas pra gerarmos a matriz G
G = GPC.G;

%% Control Solution (Analitica)
Ql = 0.8*eye(nu);
n_in = 1;
%% Simula

[y,xc] = filter(Num,Den,0); %Primeira iteracao
% x = zeros(length(H),1);
Tsim = 15;

ref = 1; %referencia de controle
w = ref*ones(n,1);
yk = y;
uk = 0;uk1 = 0;
u = [];
ne = 0;

Ub = rand(n,1); %Ubase atentar para horizonte de predição.
%Repete valores para nu < n
Ub(nu:n) = Ub(nu);
% return


Uo = zeros(nu,1);
M_inv = eye(nu*n_in);
n_ones = -1*ones(n_in*(nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;


L = [-Den(2) Num(2:end)]; % Linha do modelo SISO
Xb = zeros(n,1); %Xbase
Xo = zeros(n,1);

xk = 0;
ref = 1;
pert = 0;

for k=0:Ts:Tsim
   
    %S Model
    
    %x = xb + xo
    %y = yb + yo
    
    
    xk = [xk uk uk1]*L'; %Modelo
    
    uk1 = uk; %atualiza
    
    ne = (yk - xk);%n(t) = n(t-1) + e(t) (erro integrado)
    ni = ne*ones(n,1); %n(t+j) = n(t+j-1)   
    
    Ui = [uk; Ub(1:n-1)];
    
    xb = xk ;
    for i=1:n
        xb = [xb Ub(i) Ui(i)] *L';
        Xb(i) = xb;
    end

    
    Yb = Xb + ni;
    w = ref*ones(n,1);
    
    u0 = [uk; zeros(n_in*(nu-1),1)];
    Uc = M_inv*Ub(1:nu)-u0;
    
    %Parte desta solução pode ser computada apenas uma vez
    Uo = inv(G'*G+M_inv'*Ql*M_inv)*(G'*(w-Yb) - M_inv'*Ql*Uc);
    
    
    U = Uo + Ub(1:nu);
    uk = U(1);
    
    if(k>10)
%         pert = -0.5;
    end
    

    [yk,xc] = filter(Num(2:end),Den,uk+pert,xc); %planta
    y = [y yk]; %apenas pra armazenar e plotar
    u = [u uk];
    
    
    
end
%% 
t = 0:Ts:Tsim;
subplot(2,1,1)
plot(t,y(1:end-1),'b','linewidth',1.5); hold on;
plot(t,ref*ones(1,length(t)),'--black')
legend('y','referencia')
grid on;
title('Sistema Malha Fechada EPSAC')
subplot(2,1,2)
plot(t,u,'b','linewidth',1.5)
grid on
