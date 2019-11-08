% OBSERVACOES SOBRE A TAREFA
% otimiza_dU -> ajuste das restricoes
% gpc_tf2ss -> transforma tf para modelo GPC
% gpc_sfun -> s-function do gps. entra com matrizes do modelo GPC
%   |-> entrdadas do bloco = [leitura da planta, referencias]
clear all;close all;clc
%% MODELAGEM MIMO
Ts = 0.03; %0.03 pro mimo
z = tf('z',Ts);
s = tf('s');
H11 = tf(1,[0.7 1]);
H12 = tf(5,[0.3 1]);
H21 = tf(1,[0.5 1]);
H22 = tf(2,[0.4 1]);
Pn = [H11 H12;H21 H22];
% Pn = H11;
% Pn = H12
Pz = c2d(Pn,Ts);

% Pz = tf([0 0.4 0.6],[1 -0.8 0],Ts);
% P = Pn*exp(-0.2*s); %Insere Atraso
% P = Pn;
P = Pz
delta = [1 -1];
C = [ 1 0 0 0;1 0 0 0]';
GPC = gpc_tf2ss_U(Pz,[],delta,3,1)
% return
%% GPC - Especificacoes
lambda = 0.8;
Ql = lambda*eye(GPC.nu*GPC.nin);
Qd = 1*eye(GPC.n*GPC.nout);

Tsim = 35;
sim('simula_gpc_generic_U')

%% Plots
subplot(2,1,1)
plot(y.time,y.data)
grid on;hold on;
% legend('y_1','y_2')
% if(first == 0)
% plot(r.time,r.data,'--black')
% end
subplot(2,1,2)
plot(u.time,u.data)
% legend('u_1','u_2')
grid on;


figure(1)
subplot(2,1,1)
titel = sprintf("Planta Nominal com Restricoes | n=%d, nu=%d, \\lambda= %.2f",GPC.n,GPC.nu,lambda);
title(titel)
plot(r.time,r.data,'--black')
% legend('y_1','y_2','r_1','r_2')