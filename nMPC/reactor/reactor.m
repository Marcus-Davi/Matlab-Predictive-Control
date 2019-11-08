%Calcula uma iteracao da funcao nao linear do tanque
% h0 -> condicao inicial
% u -> vazao entrada
% Ts -> Amostragem
function x = reactor(x0,u,Ts)
Q = 0.03333;
V = 1;
k1 = 10;k2 = 10;

x = x0 + Ts*(Q/V)*(u-x0) - Ts*k1*x0/(k2*x0+1)^2;


end