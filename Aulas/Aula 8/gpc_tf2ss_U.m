%Pz -> Planta discretizada
%C -> Polinomio de robustez (n_colunas = n_saidas, n_linhas = n_estados+1)
%delta -> filtro 
%n -> horizonte de predição
%nu -> horizonte de controle
function GPC = gpc_tf2ss2(Pz,C,delta,n,nu)
[Num,Den] = tfdata(Pz);
n_sys = length(Num); %dimensao do sistema

if(isempty(C))
   noC = 1; %flag
else
    noC = 0;
end

A = [];B = []; H = []; D = [];
for i=1:n_sys %row
    d{i} = delta;
    Bk = [];
    for j=1:n_sys %coln
        d{i} = conv(d{i},Den{i,j});         
         Bi{i,j} = conv(delta,Num{i,j});%Num{i,j};%conv(delta,Num{i,j}); %B~
         for k=1:n_sys
             if(k~=j)
             Bi{i,j} = conv(Bi{i,j},Den{k,j});
             end
         end
         Bi{i,j} = [Bi{i,j}(2:end) 0]';
         Bk = [Bk Bi{i,j}];
    end
    ni = length(d{i})-1; %ordem do sistema i
    Ai{i} = zeros(ni);
    Ai{i}(:,1) = -d{i}(2:end)';
    Ai{i} = Ai{i} + diag(ones(ni-1,1),1);
    Ai{i} = blkdiag(Ai{i},0); %estado adicional
    B = [B; Bk];
    Hi{i} = zeros(1,ni+1); %estado adicional
    Hi{i}(1) = 1;    
    if(noC)
        C = zeros(ni+1,1);
        C(1) = 1;
%         dn = deconv(d{i},delta); %reproduz a dinamica da planta
%         C = conv(dn,[1 0])'; %alpha = 0% Pz = c2d(Pn,Ts);

        Di{i} = C(2:end) -d{i}(2:end)';
    else
        Di{i} = C(2:end,i) -d{i}(2:end)';
    end
    Di{i} = [Di{i} ;0]; %estado adicional
    
end

for i=1:n_sys
    A = blkdiag(A,Ai{i});
    H = blkdiag(H,Hi{i});
    D = blkdiag(D,Di{i});
end
%Matrizes GPC
G = [];F = [];E=[];
siz = length(H*B);
HB_old = [];
for i = 1:n
    HAB = H*A^(i-1)*B;
    M = [HAB HB_old repmat(zeros(siz,siz),[1 n-i])];
    G = [G;M];
    F = [F;H*A^(i)];
    E = [E;H*A^(i-1)*D];   
    HB_old = [HAB HB_old];
end

G = G(:,1:nu*n_sys);

%REVISAR CASO MIMO!! SOMAR BLOCOS 
if(nu < n)
for i=2:n
   G(i,nu) = G(i,nu) + G(i-1,nu);
end
end


% FAZER TRUNCAMENTO COM SOMATORIO

% Pré-alocação das matrizes G, F e E
% G           = zeros(N,N);
% F            = zeros(N,natil);
% E            = zeros(N,1);

% Matriz S armazena a soma dos termos h
% S            = 0;
% 
% for i=1:N
%     
%     if i==1
%         G(1,1)    = H*BB;
%     else
%         G(i,:)      = [H*(AA^(i-1))*BB G(i-1,1:end-1)];
%     end
%     F(i,:)          = H*(AA^i);
% 
%     E(i,:)          = H*(AA^(i-1))*DD;  
%     
%     S(i+1)       = S(i)+G(i,1);              
%     
% end
% 
% S          = S(2:end);
% 
% % Matriz Gnu
% Gnu                           = G(:,1:Nu);
% Gnu(Nu:end,end)     = S(1:N-Nu+1);



%Forma estrutura de dados do GPC
GPC.A = A;
GPC.B = B;
GPC.H = H;
GPC.D = D;
GPC.G = G;
GPC.F = F;
GPC.E = E;
GPC.n = n;
GPC.nu = nu;
GPC.nin = n_sys;
GPC.nout = n_sys;
GPC.Ts = Pz.Ts;

end




% M é matriz de denominadores de Pz %DIFICIL!
% function d = get_mmc(M,n)
% for i=1:n^2
%     c(i) = poly2sym(M(i,:))
% end
% 
% 
% end