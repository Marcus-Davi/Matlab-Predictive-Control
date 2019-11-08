clear;close all;clc
Ts = 1;
Tsim = 100;
time = 0:Ts:Tsim;
%% Parametros Planta
Qnl = 0.03333;
Vnl = 1;
K1nl = 10;K2nl = 10;


%% Parametros Controle
n = 5;
nu = 5;
Ql = 0.5*eye(nu);

Yb = zeros(n,1); %Xbase

yk = 0;
uk = 0;
yk1 = yk;
uk1 = uk;

n_in = 1;


M_inv = eye(nu*n_in);
n_ones = -1*ones(n_in*(nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
%% Modelo
YK = [];
UK = [];
pert = 0;
ref = 1;
for i=time 
    
        if(i == Tsim/2)
           pert = 100; 
        end
        
        if(i > Tsim/3)
           ref = 5; 
        end
        
        
        ykest = reactor(yk,uk,Ts); % RODA A PLANTA
        
         yk = reactor_uncertain(yk,uk+pert,Ts); % RODA A PLANTA
   
   e = yk-ykest; %erro de modelagem
   
    u = uk*ones(nu,1);
%     y = cost(ref,yk,u,Ts,n);
   w = ref*ones(n,1);
   
   
   
   uk_newton = uk;
   for iter=1:10
   
   yb = yk;
   x = zeros(n,1);
    for j=1:n
     yb = reactor(yb,uk_newton,Ts) + e;   %incluir erro
     x(j) = yb;
    end
      
   jac = (Ts*Qnl/Vnl)*ones(n,1);
   residuals = Ql*(x-w);
   increment = -inv(jac'*jac)*jac'*residuals ;
%     increment = -1/jac*residuals;

    uk_newton = uk_newton + increment;
   
   end
   
   
   
%     X = fmincon(@(u) cost(ref,yk,e,u,Ts,n),u,[],[]);
    
    
%     uk = X(1);
    
    uk = uk_newton
    
    if(uk > 100)
        uk = 100;
    end
%     uk = 0.1392
   YK = [YK yk];
   UK = [UK uk];
   i
end

subplot(2,1,1)
plot(time,YK)
grid on
subplot(2,1,2)
plot(time,UK)
grid on;
xlabel('Tempo [s]')
ylabel('C [mol]')


%% Cost

function [y,naoimporta] = cost(ref,yk,e,u,Ts,n)
x = zeros(n,1);

yb = yk;
    for i=1:n
     yb = reactor(yb,u(i),Ts) + e;   %incluir erro
     x(i) = yb;
    end
    
    y = norm(ref-x)^2 + 0.9*norm(u)^2;
    
naoimporta =[];
end

function r = residual(m,b,x,y)
    r =  m*exp(x)+b - y;
end

function j = jacobian(m,b,x)
   j =  [exp(x);1];
end
