clear;close all;clc
%% Planta
Ts = 1;
Num = [1 -0.8 0];
Den = [0 0.4 0.6] ;

delta = [1 -1];
nB = length(Den)-1;
nA = length(Num)-1;
Atil = conv(Num,delta);

A = [-Atil(2) 1;-Atil(3) 0];
B = [Den(2) Den(3)]';
P = tf(Den,Num,Ts)

LEG = {};
%% Matrizes do sistema
for alpha = 0:0.1:0.5

c1 = [1 -alpha];
c2 = [1 -alpha];

C = conv(c1,c2)
D = [C(2) - Atil(2);C(3) - Atil(3)];
H = [1 0];

%Horizonte
n = 3;
nstate = length(A);

G = [];F = [];E=[];
siz = length(H*B);
HB_old = [];
for i = 1:n
    HAB = H*A^(i-1)*B;
    M = [HAB HB_old repmat(zeros(siz,siz),[1 n-i])];
    G = [G;M];
    F = [F;H*A^(i)];
    E = [E;H*A^(i-1)*D];   
    HB_old = [HAB HB_old];
end
G
F
E
%% Solution
nu = 3;
G = G(:,1:nu);
Ql = 1*eye(nu);
K = inv(G'*G + Ql)*G';
%dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)
Kr = K(1,:) * ones(n,1);
Kf = K(1,:) * F;
Ke = K(1,:) * E;
%% Simula
Tsim = 20;
sim('simula_ss')


%% Plots
figure(1)
subplot(2,1,1)
plot(y.time,y.data)
grid on;hold on;
subplot(2,1,2)
plot(u.time,u.data)
grid on;hold on;
concat = sprintf('alpha = %.2f',alpha);
LEG = [LEG concat];
end

%% Legends
figure(1)
subplot(2,1,1)
legend(LEG)

