clear all;close all;clc
addpath('../GPC_LIB')
%% MODELAGEM MIMO
Ts = 0.03;
z = tf('z',Ts);
s = tf('s');
H11 = tf(1,[0.7 1]);
H12 = tf(5,[0.3 1]);
H21 = tf(1,[0.5 1]);
H22 = tf(2,[0.4 1]);
Pn = [H11 H12];
% Pn = [H11 H12;H21 H22];
% Pn = H11
Pz = c2d(Pn,Ts);
% P = Pn*exp(-0.2*s); %Insere Atraso
P = Pn;
delta = [1 -1];
alpha1 = 0;
alpha2 = 0;
c1 = [1 -alpha1];
c2 = [1 -alpha2];


C = [ 1 0 0 0;1 0 0 0]';

% [A,B,H,D] = gpc_tf2ss_U(Pz,C,delta)

n = 3; %Horizonte AJUSTAR AQUI
nu = 3; %Horizonte controle AJUSTAR AQUI


GPC = gpc_tf2ss_midRange(Pz,C,delta,n,nu);

%% Solução Analítica
Qdelta = eye(n*GPC.nout);
lambda = 0.8;
Q1 = lambda*eye(GPC.nu);
Q2 = lambda*eye(GPC.nu);
Ur = 0.1*ones(n,1);
w = 1*ones(n,1);
G1 = GPC.G1;
G2 = GPC.G2;
H = [G1'*Qdelta*G1+Q1 G1'*Qdelta*G2;G2'*Qdelta*G1 G2'*Qdelta*G2+Q2];
b = [G1'*Qdelta G2'*Qdelta];
K = inv(H)*b';

nin = GPC.nin;
Kf = K(nin,:)*GPC.F;
Ke = K(nin,:)*GPC.E;
Kr = K(nin,:)*w;
K2 = K(nin,:,1)





