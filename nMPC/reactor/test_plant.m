clear all;close all;clc
Ts = 0.1;
Tsim = 500;
time = 0:Ts:Tsim;

Qnl = 0.03333;
Vnl = 1;
K1nl = 10;K2nl = 10;

x = 0;
u = 8;
X = [];
for i=time
    x = reactor(x,u,Ts);
    
    
    
    
   
    X = [X x];
end

plot(time,X)