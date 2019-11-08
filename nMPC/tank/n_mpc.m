clear;close all;clc
Ts = 60;
Tsim = Ts*100;
time = 0:Ts:Tsim;

%% Parametros Controle
n = 4;
nu = n;
Ql = 0*eye(nu);

Yb = zeros(n,1); %Xbase


n_in = 1;


% M_inv = eye(nu*n_in);
% n_ones = -1*ones(n_in*(nu-1),1);
% n_diag = diag(n_ones,-n_in);
% M_inv = M_inv+n_diag;
%% Modelo
C = [1 -0.7]; %filtro do erro
ef = 0; %erro filtrado
YK = [];
UK = [];
NOISE = [];
pert = 0;



yk = 1; %inicial
uk = 0.1392; %inicial

ykest = yk;
ref = 1;

for i=time 
    time_s = i/60;
        if(time_s  > 5 && time_s < 15 )
        pert = -0.1;
        elseif(time_s > 20)
           ref = 2; 
        end
    
%         ykest = tank(ykest,uk,Ts); % PARAELO
      ykest = tank(yk,uk,Ts); % SERIE-PARALELO
      
        noise = 0.4*randn;
        
         yk = tank_incerto(yk,uk+pert,Ts); % RODA A PLANTA
   
        e = yk-ykest + noise; %medição!
        ef = e - C(2)*ef; %erro filtrado por C
        
       
        
  
        %TODO GENERALIZAR USANDO DECONV!
        F1 = sum(C);
        F2 = F1;
        F3 = F1;
        
        nf(1) = F1*ef;
        nf(2) = F2*ef;
        nf(3) = F3*ef;
        
        
        
   
         u = uk*ones(nu,1);
         w = ref*ones(n,1); 
   
%    uk_newton = uk*ones(n,1);
%    for iter=1:5
%    
% %    yb = yk; %SERIE PARALELO
% yb = ykest; %PARALELO
%    x = zeros(n,1);
%     for j=1:n
%      yb = tank(yb,uk_newton(j),Ts) + e;   %incluir erro
%      x(j) = yb;
%     end
%       
% %    jac = (Ts/10)*ones(n,1);
% jac = (Ts/10)*eye(n);
%    residuals = (x-w);
%    increment = -inv(jac'*jac)*jac'*residuals ;
% %     increment = -1/jac*residuals;
% 
%     uk_newton = uk_newton + increment;
%    
%    end


    nf = e*ones(n,1);
    
%     X = fmincon(@(u) cost(ref,ykest,nf,u,uk,Ts,n,nu,C,Ql),u,[],[]); %PARALELO
      X = fmincon(@(u) cost(ref,yk,nf,u,uk,Ts,n,nu,C,Ql),u,[],[]); %SERIE-PARALELO
   
     uk = X(1);
%     uk = uk_newton(1);
%     uk = 0.1392
   YK = [YK yk];
   UK = [UK uk];
   NOISE = [NOISE noise];
end
time = time/60;
plot(time,YK)
hold on
plot(time,UK)
grid on;
xlabel('Tempo [m]')
ylabel('Altura Tanque [m]')
% figure
% plot(time,NOISE)
% grid on

%% Cost

function [y,naoimporta] = cost(ref,yk,nf,u,uk,Ts,n,nu,C,Ql)
x = zeros(n,1);

[~,n_in] = size(uk);

M_inv = eye(nu*n_in);
n_ones = -1*ones(n_in*(nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;

nlin = length(M_inv);
u0 = [uk; zeros(nlin-1,1)];


yb = yk;
    for i=1:n
     yb = tank(yb,u(i),Ts) + nf(i);   %incluir erro
     x(i) = yb;
    end
    
    y = norm(ref-x)^2 + (M_inv*u-u0)'*Ql*(M_inv*u-u0);
    
naoimporta =[];
end

