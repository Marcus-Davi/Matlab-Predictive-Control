%Computa G variando a amostragem 
function G = get_G_var(IC,model,increment,e,n,nu,Ts)
x0 = IC.x0;
uin = IC.u0;

nin = length(uin);
nout = length(x0);

% USAR BLOCOS!
%blkzise = nout x nin
blksize = [nout nin];
H = zeros(n*nout,nin);
DH = H;
u0 = repmat(uin,1,n);
    for i=1:nin %entradas
        h0 = x0;
        dh0 = x0;
        
     u0_inc = u0;
     u0_inc(i,1) = u0(i,1) + increment;

     for j=1:n %predição
         
    
    h = model(h0,u0(:,j),j*Ts) + e(:,j); 
    dh = model(dh0,u0_inc(:,j),j*Ts) + e(:,j);  
    
    h0 = h;
    dh0 = dh;
   
    H = set_block(H,j,i,[nout,1],h);
    DH = set_block(DH,j,i,[nout,1],dh);
    
     end
   end
%aqui h é o primeiro bloco de G
delta = (DH-H)/increment;
deltaH = delta;
% for i=2:n
%     blk = get_block(delta,i,1,blksize) - get_block(delta,i-1,1,blksize);
%     deltaH = set_block(deltaH,i,1,blksize,blk);
% %    HH(i) = H(i)-H(i-1);  %derivada
% end

G = [];
HB_old = [];
for i = 1:n
    HAB = get_block(deltaH,i,1,blksize);
    M = [HAB HB_old repmat(zeros(blksize),[1 n-i])];
    G = [G;M]; 
    HABi{i} = HAB;
    HB_old = [HAB HB_old];
end

for i=2:n
   HABi{i} = HABi{i} + HABi{i-1} ;
end

G = G(:,1:nu*nin);

%REVISAR CASO MIMO!! SOMAR BLOCOS.
if(nu < n)
       Gu = repmat(zeros(nout,nin),[nu-1 1]);
for i=1:(n-nu+1)
   Gu = [Gu;HABi{i}];
end
G(:,(nu*nin-(nin-1)):end) = Gu;
end



end



