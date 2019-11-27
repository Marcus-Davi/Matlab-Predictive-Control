clear;clc;close all;
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