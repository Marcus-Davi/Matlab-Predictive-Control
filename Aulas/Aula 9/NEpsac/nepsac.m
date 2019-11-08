clear;close all;clc
Ts = 60;
Tsim = Ts*50;
time = 0:Ts:Tsim;


%Iterações

%% Parametros Controle
n = 3;
nu = 3;
Ql = 0.8*eye(nu);

Uo = zeros(nu,1);
% Xb = zeros(n,1); %Xbase
Yb = zeros(n,1); %Xbase

yk = 1;
uk = 0.1392;

n_in = 1;
ref = 50;

M_inv = eye(nu*n_in);
n_ones = -1*ones(n_in*(nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;

%0 -> PARAELO | 1-> SERIE_PARALELO    <----- ESCOLHA AQUI ###
SERIE_PARALELO = 1; 
%% Modelo
YK = [];
UK = [];
pert = 0;
xk = yk;
for i=time 
   
    if(i  > 500)
        pert = 0.5;
    end
    
    if(SERIE_PARALELO)
        xk = tank(yk,uk,Ts); %SERIE PARALELO          
    else
     xk = tank(xk,uk,Ts); %PARALELO (CUIDADO COM RAIZ NEGATIVA P ESSE SISTEMA!)
    end
    
     yk = tank(yk,uk+pert,Ts); % RODA A PLANTA
     
     
     
    ne = (yk - xk);%n(t) = y(t) - x(t)
    ni = ne*ones(n,1); %n(t+j) = n(t+j-1) 
    
    Ub = uk*ones(n,1); %Ubase
    Ub(nu:n) = Ub(nu);
    
    
%     yb = xk; % PARALELO
if(SERIE_PARALELO)
    yb = yk; % SERIE-PARALELO
    for k=1:n
     yb = tank(yb,Ub(k),Ts) + ni(k);  
     Yb(k) = yb;
    end
    
    else
        yb = xk;
    for k=1:n
        yb = tank(yb,Ub(k),Ts) ;
        Yb(k) = yb;
    end
    Yb = Yb + ne;
    
end
    
    IC.x0 = yk;
    IC.u0 = Ub(1);
    
   G = get_G(IC,@tank,0.00001,n,nu,Ts);
   

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

