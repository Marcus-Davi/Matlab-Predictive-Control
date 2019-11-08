clear;close all;clc
s = tf('s');
Ts = 0.1;
QSI = 1;
Pn = tf(1,[1 2*QSI 1]);
Pz = c2d(Pn,Ts);
% P = Pn*exp(-0.2*s);
P = Pn;
[Num,Den] = tfdata(Pz,'v');
% Num = [1 -0.8 0];
% Den = [0 0.4 0.6] ;

delta = [1 -1];
nB = length(Num)-1;
nA = length(Den)-1;
Atil = conv(Den,delta)';

A = [-Atil(2:end) [1 0 0]' [0 1 0]'];
B = [Num(2:end) 0]';

LEG = {};
%% Matrizes do sistema
alpha = 0.;
c = [1 -alpha];
% 

C = conv(conv(c,c),c)';
% C(2) = -0.5;
D = C(2:end) - Atil(2:end);
H = [1 0 0];

%Horizonte

% n = 10;
% for n = 10:1:20
for n = 1
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
G;
F;
E;
% return
%% Solution
nu = n;
G = G(:,1:nu);
lambda = 0.0;
Ql = lambda*eye(nu);
K = inv(G'*G + Ql)*G';
%dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)
Kr = K(1,:) * ones(n,1);
Kf = K(1,:) * F;
Ke = K(1,:) * E;
%% Simula
Tsim = 20;
sim('simula_ajuste_siso')


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

%% Legends
figure(1)
subplot(2,1,1)
% ylim([-10 10])
% legend(LEG)
pause(0.5)
n
end
