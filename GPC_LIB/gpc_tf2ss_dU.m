%Pz -> Planta discretizada
%C -> Polinomio de robustez (n_colunas = n_saidas, n_linhas = n_estados+1)
%delta -> filtro 
%n -> horizonte de predição
%nu -> horizonte de controle
function GPC = gpc_tf2ss_dU(Pz,C,delta,n,nu)
[Num,Den] = tfdata(Pz);
[nout,nin] = size(Num); %dimensao do sistema

if(isempty(C))
   noC = 1; %flag
else
    noC = 0;
end

A = [];B = []; H = []; D = [];
for i=1:nout %row
    d{i} = delta;
    Bk = [];
    for j=1:nin %coln
        d{i} = conv(d{i},Den{i,j});         
         Bi{i,j} = Num{i,j};
         for k=1:nin
             if(k~=j)
             Bi{i,j} = conv(Bi{i,j},Den{i,k});
             end
         end
         Bi{i,j} = [Bi{i,j}(2:end) 0]';
         Bk = [Bk Bi{i,j}];
    end
    ni = length(d{i})-1; %ordem do sistema i
    Ai{i} = zeros(ni);
    Ai{i}(:,1) = -d{i}(2:end)';
    Ai{i} = Ai{i} + diag(ones(ni-1,1),1);
    B = [B; Bk];
    Hi{i} = zeros(1,ni);
    Hi{i}(1) = 1;    
    if(noC)
        C = zeros(ni+1,1);
        C(1) = 1;
%         dn = deconv(d{i},delta); %reproduz a dinamica da planta
%         C = conv(dn,[1 0])'; %alpha = 0
        Di{i} = C(2:end) -d{i}(2:end)';
    else
        Di{i} = C(2:end,i) -d{i}(2:end)';
        
    end
    
end

[~,sA] = size(Ai);
[~,sH] = size(Hi);

for i=1:sA
   A = blkdiag(A,Ai{i}); 
end

for i=1:sH
    H = blkdiag(H,Hi{i});
    D = blkdiag(D,Di{i});  
end
%Matrizes GPC
G = [];F = [];E=[];
HB_old = [];
for i = 1:n
    HAB = H*A^(i-1)*B;
    M = [HAB HB_old repmat(zeros(nout,nin),[1 n-i])];
    G = [G;M];
    F = [F;H*A^(i)];
    E = [E;H*A^(i-1)*D];   
    HB_old = [HAB HB_old];
end

G = G(:,1:nu*nin);

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
GPC.nin = nin;
GPC.nout = nout;
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