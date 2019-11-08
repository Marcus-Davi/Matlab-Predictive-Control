clear;close all;clc
A = [1 -0.8];
B = [0.4 0.6];
C = 1;
delta = [1 -1];
nB = length(B);
nA = length(A);

Atil = conv(A,delta);
A1 = -Atil(2:end);
B1 = B(2);
%AY(t) = BU(t) + C/delta E(t)
n =4  ;
%Y(1) = ((B(2:end,:)*p)'*U(2:end))* (A(1,:)*p);
% Y(t) = 1.8 y(t-1) - 0.8 y(t-2) + 0.4 u(t-1) + 0.6 u(t-2)
% Y(t+1) = 1.8 y(t) - 0.8 y(t) + 0.4 u(t-1) + 0.6 u(t-2)
M = zeros(n,n); %delta [u(t) u(t+1) +  ...]
M1 = zeros(n,1); %delta u(t-1)
M2 = zeros(n,2); % y(t)

for i=1:n
    if(i==1)
        m = B(1);
        An = A1;
        Bn = B1;
        M(i,i) = m;
    elseif(i==nA) %varia com a ordem?
       m = A1(1) * m + B1;
    M = add_and_shift(M,m,i);
    An = [A1(1)*An(1) + An(2) A1(1)*An(2)]; %regra 1
    Bn = A1(1)*Bn;
    else
        m = A1(1) * m + B(1)*A1(2);
        M = add_and_shift(M,m,i);
   An = [A1(1)*An(1) + An(2) A1(1)*An(2) + A1(2)^2]; %regra geral (melhorar?)
   Bn = A1(1)*Bn+(B(2)*A1(2))
    end
    
    
    M2(i,:) = An;
    M1(i,:) = Bn;
end

M
M1
M2