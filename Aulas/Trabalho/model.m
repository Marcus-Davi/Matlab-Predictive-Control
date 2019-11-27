clear;clc;close all;
addpath('../GPC_Funcs_dU');
%% System Model
Ts = 0.1;
Tau =  0.2;
s = tf('s');
z = tf('z',Ts);
H = tf(1,[Tau 1]);
% H = tf(1,[1 0.05]);
L = 1;
P = H*exp(-L*s);

Pz = c2d(H,Ts);
d = round(L/Ts);
C = [1 0 0]';
delta = [1 -1];

%% Control Solution
n = 3; nu = 3;
GPC = gpc_tf2ss_dU(Pz,C,delta,n,nu);
G = GPC.G;
A = GPC.A;
B =GPC.B;
D = GPC.D;
H = GPC.H;

Ql = 0.8*eye(nu);
K = inv(G'*G + Ql)*G';
Kr = K(1,:) * ones(n,1);
Kf = K(1,:) * GPC.F;
Ke = K(1,:) * GPC.E;
Ad = GPC.A^d;
AdD = GPC.A^(d-1)*GPC.D;
sum = 0;
AB = [0;0];
for k=1:d
   AB = [AB GPC.A^(k-1)*GPC.B] %AB -> Filtro FIR
end 

%% Simula
Tsim = 10;
sim('simula_model')
subplot(2,1,1)
plot(y.time,y.data,'Linewidth',2)
legend('Saída')
grid on;
subplot(2,1,2)
plot(u.time,u.data,'r','Linewidth',2)
legend('Sinal de Controle')
grid on;


