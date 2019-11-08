clear;close all;clc
Ts = 1;
A = [1 -0.8 0];
B = [0 0.4 0.6] ;
C = 1;
% d = 0;


% z = tf('z',Ts);
% B = (0.4*z^-1 + 0.6*z^-2)*z^-d;
% A = 1 - 0.8*z^-1 + 0.3*z^-2;
% P = minreal(B/A);
% [B,A] = tfdata(P,'v')

delta = [1 -1];
nB = length(B)-1;
nA = length(A)-1;
Atil = conv(A,delta);

nB = 1
P = minreal(tf(B,A,Ts))


% Calcula divisao polinomial
n = 3;
Bi = [1 zeros(1,nA+n)];
[E,Q] = deconv(Bi,Atil);
%Resto sempre tera 2 termos (os 2 ultimos)
EB = conv(B,E);
G = zeros(n,n);
F = zeros(n,nA+1);
Gp = zeros(n,1);

for i=1:n
   G(i,:) = EB(2:n+1);
   b = [1 zeros(1,i+nA)];
   [q,f] = deconv(b,Atil);
   F(i,:) = f(end-(nA):end);
   eb = conv(q,B);
   Gp(i,1) = eb(find(eb,1,'last'));
end
G = ajeita_G(G)
F
Gp
% return
%% A partir daqui temos G, Gp, F
% Lei de controle
nu = 3;
G = G(:,1:nu);
Ql = 1*eye(nu);
K = inv(G'*G + Ql)*G';
%dU = -(G'G + Ql)G' (w - f) = K (w - f) ; f = FY(t) + GpU(t-1)
% como nosso w = [1 ... 1]' * r(t), dU = Kr*r - Kfy(t) 0 KgpdU(t-1)
Kr = K(1,:) * ones(n,1);
Kf = K(1,:) * F;
Kgp = K(1,:) * Gp;
% return
sim('Simula')
%% Plota
time = u.time;
subplot(2,1,1)
plot(time,y.data)
grid on;
subplot(2,1,2)
plot(time,u.data)
grid on;