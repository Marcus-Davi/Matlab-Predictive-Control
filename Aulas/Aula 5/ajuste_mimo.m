clear;close all;clc
Ts = 0.03;
z = tf('z',Ts);
s = tf('s');
H11 = tf(1,[0.7 1]);
H12 = tf(5,[0.3 1]);
H21 = tf(1,[0.5 1]);
H22 = tf(2,[0.4 1]);
Pn = [H11 H12;H21 H22];
Pz = c2d(Pn,Ts);

P = Pn*exp(-0.2*s); %Insere Atraso
% P = Pn;
delta = [1 -1];
d11 = Pz(1,1).den{1};n11 = Pz(1,1).num{1};
d12 = Pz(1,2).den{1};n12 = Pz(1,2).num{1};
% d13 = Pz(1,3).den{1};n13 = Pz(1,3).num{1};

d21 = Pz(2,1).den{1};n21 = Pz(2,1).num{1};
d22 = Pz(2,2).den{1};n22 = Pz(2,2).num{1};
% d23 = Pz(2,3).den{1};n23 = Pz(2,3).num{1};


% d31 = Pz(3,1).den{1};n31 = Pz(3,1).num{1};
% d32 = Pz(3,2).den{1};n32 = Pz(3,2).num{1};
% d33 = Pz(3,3).den{1};n33 = Pz(3,3).num{1};

A1til = conv(conv(delta,d11),d12)';
A2til = conv(conv(delta,d21),d22)';

B11 = conv(n11,d12)';
B12 = conv(n12,d11)';
B21 = conv(n21,d22)';
B22 = conv(n22,d21)';


A1 = [-A1til(2:end) [1 0 0]' [0 1 0]'];
B1 = [B11 B12];
B1 = [B1(2:end,:);0 0]; %ajeita

A2 = [-A2til(2:end) [1 0 0]' [0 1 0]'];
B2 = [B21 B22];
B2 = [B2(2:end,:);0 0]; %ajeita

H1 = [1 0 0];
H2 = [1 0 0];

nA = length(A1);
A = [A1 zeros(nA);zeros(nA) A2];
B = [B1;B2];
H = [H1 zeros(1,3);zeros(1,3) H2];
% LEG = {};
first = 0;
%for alpha = [0]
for alpha = 0.9:0.01:0.92
c = [1 -alpha];
% C1 = [1 0 0 0]';
C1 = conv(conv(c,c),c)';
D1 = C1(2:end) - A1til(2:end);

% C2 = [1 0 0 0]';
C2 = C1;
D2 = C2(2:end) - A2til(2:end);
nD = length(D1);
D = [D1 zeros(nD,1);zeros(nD,1) D2];

%% GPC
n = 5; %Horizonte AJUSTAR AQUI
nu = 5; %Horizonte controle AJUSTAR AQUI
lambda = 0.01; %AJUSTAR AQUI
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

[~,nin] = size(B); %n entradas
G = G(:,1:nu*nin);
Ql = lambda*eye(nu*nin);
K = inv(G'*G + Ql)*G';
%dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)
Ku = K(1:nin,:);
Kr = 0;
for i=1:n
Kr = Kr + Ku(1:nin,nin*i-1:nin*i); 
end
Kf = Ku * F;
Ke = Ku * E;


Tsim = 8;
sim('simula_ajuste_mimo')

%% Plots

subplot(2,1,1)
plot(y.time,y.data)
grid on; hold on;
% if(first == 0)
% plot(r.time,r.data,'--black')
% end
subplot(2,1,2)
plot(u.time,u.data)
grid on;hold on;
end
figure(1)
subplot(2,1,1)
titel = sprintf("Planta com Atraso | n=%d, nu=%d, \\lambda= %.2f",n,nu,lambda);
title(titel)
plot(r.time,r.data,'--black')
% legend(LEG)