clear;close all;clc
Ts = 1;
z = tf('z',Ts);

H11 = tf(0.0420,[1 -0.9580],Ts);
H12 = tf(0.4758,[1 -0.9048],Ts);
H21 = tf(0.0582,[1 -0.9418],Ts);
H22 = tf(0.1445,[1 -0.9277],Ts);
Pn = [H11 H12;H21 H22];
P = 1.2*Pn; 
delta = [1 -1];

d11_d12 = conv(delta,H11.den{1});
d11_d12 = conv(d11_d12,H12.den{1}); %A1`

d21_d22 = conv(delta,H21.den{1});
d21_d22 = conv(d21_d22,H22.den{1}); %A2`

n11_d12 = conv(H11.num{1},H12.den{1});%B11
n12_d11 = conv(H12.num{1},H11.den{1}); %B12

n21_d22 = conv(H21.num{1},H22.den{1}); %B21
n22_d21 = conv(H22.num{1},H21.den{1}); %B22

A1til = d11_d12;
A2til = d21_d22;

B11 = n11_d12;
B12 = n12_d11;
B21 = n21_d22;
B22 = n22_d21;



A1 = [-A1til(2) 1 0;-A1til(3) 0 1;-A1til(4) 0 0];
B1 = [B11' B12'];
B1 = [B1(2:end,:);0 0]; %ajeita


LEG = {};
first = 0;
for alpha = 0
c = [1 -alpha];
% C1 = [1 0 0 0]';
C1 = conv(conv(c,c),c)';

D1 = C1(2:end) - A1til(2:end)';

A2 = [-A2til(2) 1 0;-A2til(3) 0 1;-A2til(4) 0 0];
B2 = [B21' B22'];
B2 = [B2(2:end,:);0 0]; %ajeita
% C2 = [1 0 0 0]';
C2 = C1;
D2 = C2(2:end) - A2til(2:end)';


H1 = [1 0 0];
H2 = [1 0 0];

nA = length(A1);
nD = length(D1);
A = [A1 zeros(nA);zeros(nA) A2];
B = [B1;B2];
H = [H1 zeros(1,3);zeros(1,3) H2];
D = [D1 zeros(nD,1);zeros(nD,1) D2];

%% GPC
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
nu = n;
[~,nin] = size(B); %n entradas
G = G(:,1:nu*nin);
Ql = 0.05*eye(nu*nin);
K = inv(G'*G + Ql)*G';
%dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)
Ku = K(1:nin,:);
Kr = Ku(1:2,1:2) + Ku(1:2,3:4) + Ku(1:2,5:6) %+ Ku(1:2,7:8) + Ku(1:2,9:10); %Automatizar
Kf = Ku * F;
Ke = Ku * E;


Tsim = 100;
sim('simula_mimo')

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
concat = sprintf('alpha = %.2f',alpha);
LEG = [LEG concat];
end

figure(1)
subplot(2,1,1)
title('Planta com 20% de incerteza no ganho')
plot(r.time,r.data,'--black')
legend(LEG)