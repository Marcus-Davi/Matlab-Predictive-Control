clear;close all;clc
%% Planta
Ts = 1;
Den = [1 -0.8 0];
Num = [0 0.4 0.6] ;

delta = [1 -1];
nB = length(Num)-1;
nA = length(Den)-1;
Atil = conv(Den,delta);

A = [-Atil(2) 1;-Atil(3) 0];
B = [Num(2) Num(3)]';
P = tf(Num,Den,Ts)

%% Projeto Controle
alpha = 0.8;
c1 = [1 -alpha];
c2 = [1 0];

C = conv(c1,c2)
D = [C(2) - Atil(2);C(3) - Atil(3)];
H = [1 0];

%Horizonte
n = 3;
nu = 3; %horizonte de controle
nstate = length(A);

G = [];F = [];E=[];
s = length(H*B);
HB_old = [];
for i = 1:n
    HAB = H*A^(i-1)*B;
    M = [HAB HB_old repmat(zeros(s,s),[1 n-i])];
    G = [G;M];
   
    F = [F;H*A^(i)];
    E = [E;H*A^(i-1)*D];
    
    HB_old = [HAB HB_old];
    
end
G
F
E


%% Control Solution (Analitica)

G = G(:,1:nu);
Ql =0.8*eye(nu);
K = inv(G'*G + Ql)*G';
%dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)
Kr = K(1,:) * ones(n,1);
Kf = K(1,:) * F;
Ke = K(1,:) * E;

%% Restricoes -> Alguns valores podem causar infactibilidade
M = tril(ones(nu));
Umax = 0.4;
Umin = -0.4;
dUmax = 0.5;
dUmin = -0.5;
Ymax = 1.5;
Ymin = -0.5;
options = optimoptions('quadprog','Display','off')
%% Simula

[y,xc] = filter(Num,Den,0); %Primeira iteracao
x = zeros(length(H),1);
Tsim = 30;

ref = 1; %referencia de controle
w = ref*ones(n,1);
yk = 0;
uk = 0;
du = [];
u = [];
pert = 0;

for k=0:Ts:Tsim
   
    e = yk - H*x;
    f = F*x + E*e;
    
    dU = K*(w-f); %solucao analitica
%     dU = quadprog(G'*G+Ql,G'*(f-w),[],[]) %sem restricoes
    
    u0 = [uk; zeros(nu-1,1)]; %apenas a uk é importante
    %UMAX
    A_con1 = [M;-M];
    B_con1 = [Umax-u0;u0-Umin];
    %DUMAX    
    A_con2 = [eye(nu);-eye(nu)];
    B_con2 = [dUmax*ones(nu,1);-dUmin*ones(nu,1)];
    %YMAX
    A_con3 = [G;-G];
    B_con3 = [Ymax-f;f-Ymin];
    
    %Aqui escolhemos quais das restrições anteriores são ustadas. Basta
    %concatenar;
    A_con = [A_con1;A_con2;A_con3];
    B_con = [B_con1;B_con2;B_con3];
    
    dU = quadprog(G'*G+Ql,G'*(f-w),[],[],[],[],[],[],[],options); %com restricoes
    
    
    % Usando Yalmip (Alternativa)
%     Y_Con = A_con*x_sdp - B_con <= 0;
%     Y_H = G'*G+Ql;
%     Y_b = G'*(f-w);
%     Y_Min = 0.5*x_sdp'*Y_H*x_sdp + Y_b'*x_sdp;
%     optimize(Y_Con,Y_Min);
%     dU = value(x_sdp);
    
    uk = dU(1) + uk;
    x = B*dU(1) + D*e + A*x;
    
    if(k > 10)
       pert = -0.1;
    end
   
    [yk,xc] = filter(Num(2:end),Den,uk+pert,xc); %planta
    y = [y yk]; %apenas pra armazenar e plotar
    u = [u uk];
    du = [du dU(1)];
    
end

t = 0:Ts:Tsim;
subplot(2,1,1)
plot(t,y(1:end-1),'b','linewidth',1.5); hold on;
plot(t,ref*ones(1,length(t)),'--black')
plot(t,Ymax*ones(1,length(t)),'--b')
plot(t,Ymin*ones(1,length(t)),'--b')
legend('y','referencia')
grid on;
title('Sistema de Controle com Restricoes')
subplot(2,1,2)
plot(t,u,'b','linewidth',1.5)
hold on;
plot(t,du,'r','linewidth',1.5);
plot(t,Umax*ones(1,length(t)),'b--')
plot(t,Umin*ones(1,length(t)),'b--')

plot(t,dUmax*ones(1,length(t)),'r--')
plot(t,dUmin*ones(1,length(t)),'r--')
legend('u','\Deltau')
grid on;