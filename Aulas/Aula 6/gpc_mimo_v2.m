% OBSERVACOES SOBRE A TAREFA
% otimiza_dU -> ajuste das restricoes
% gpc_tf2ss -> transforma tf para modelo GPC
% gpc_sfun -> s-function do gps. entra com matrizes do modelo GPC
%   |-> entrdadas do bloco = [leitura da planta, referencias]


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
Pn = [H11 H12;H21 H22];
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

n = 5; %Horizonte AJUSTAR AQUI
nu = 1; %Horizonte controle AJUSTAR AQUI


GPC = gpc_tf2ss_dU(Pz,C,delta,n,nu)
A = GPC.A;
B = GPC.B;
H = GPC.H;
D = GPC.D;



%% GPC - Especificacoes
[~,n_in ] = size(B);
lambda = 1;
Ql = lambda*eye(nu*n_in);
Ql(1,1) = 0;
Ql(2,2) = 0;
Qd = 1*eye(n*n_in);

Tsim = 3;
sim('simula_mimo_v2')

%% Plots
subplot(3,1,1)
plot(y.time,y.data)
grid on;hold on;
legend('y_1','y_2')
% if(first == 0)
% plot(r.time,r.data,'--black')
% end
subplot(3,1,2)
plot(u.time,u.data)
legend('u_1','u_2')
grid on;
subplot(3,1,3)
plot(du.time,du.data)
legend('\delta u_1','\delta u_2')
grid on;


figure(1)
subplot(3,1,1)
titel = sprintf("Planta Nominal com Restricoes | n=%d, nu=%d, \\lambda= %.2f",n,nu,lambda);
title(titel)
plot(r.time,r.data,'--black')
% legend(LEG)