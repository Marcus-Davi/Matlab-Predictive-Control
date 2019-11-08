clear;close all;clc
Ts = 60;
Tsim = 60*15;
time = 0:Ts:Tsim;


%Iterações

%% Parametros Controle
n = 3;
nu = 3;
Ql = 1*eye(nu);

Uo = zeros(nu,1);
% Xb = zeros(n,1); %Xbase
Yb = zeros(n,1); %Xbase

yk = 1;
uk = 0.1392;
yk1 = yk;
uk1 = uk;

n_in = 1;
ref = 5;


M_inv = eye(nu*n_in);
n_ones = -1*ones(n_in*(nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
%% Modelo
YK = [];
UK = [];
pert = 0;

for i=time 
   
    if(i  > 500)
        pert = 1;
    end
    
   yk = tank_incerto(yk,uk+pert,Ts); % RODA A PLANTA

   xk = tank(yk1,uk,Ts); %SERIE-PARALELO x(y) = f(y(t-1),u(t-1)
   
   yk1 = yk;
   uk1 = uk;
   
    ne = (yk - xk);%n(t) = y(t) - x(t)
    ni = ne*ones(n,1); %n(t+j) = n(t+j-1) 
    
    Ub = uk1*ones(n,1); %Ubase
    Ub(nu:n) = Ub(nu);
    
    yb = yk;
    for i=1:n
     yb = tank(yb,Ub(i),Ts) + ni(i);  
     Yb(i) = yb;
    end
    
    IC.x0 = yk;
    IC.u0 = uk;
    
%    G = get_G_old(IC,@tank,0.01,n,nu,Ts);
   G = get_G(IC,@tank,0.01,n,nu,Ts);
   K0 = inv(G'*G+M_inv'*Ql*M_inv);
   w = ref*ones(n,1);
   
    u0 = [uk; zeros(n_in*(nu-1),1)];
    Uc = M_inv*Ub(1:nu)-u0;
    
    Uo = K0*(G'*(w-Yb) - M_inv'*Ql*Uc);
    
    
    U = Uo + Ub(1:nu);
    uk = U(1);
   
   
   YK = [YK yk];
   UK = [UK uk];
end
time = time/60;
plot(time,YK)
hold on
plot(time,UK)
grid on;
xlabel('Tempo [m]')
ylabel('Altura Tanque [m]')

