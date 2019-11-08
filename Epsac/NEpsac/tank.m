%Calcula uma iteracao da funcao nao linear do tanque
% h0 -> condicao inicial
% u -> vazao entrada
% Ts -> Amostragem
function h = tank(h0,u,Ts)

g = 9.81;
A = 10;
a = 0.01*pi;
k1 = sqrt(2*g)*a;

h = h0 + Ts*u/A - k1*sqrt(h0)*Ts/A;


end