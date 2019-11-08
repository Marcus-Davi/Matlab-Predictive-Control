function U = otimiza_U(G,Ql,Qd,f,w,nu,n,n_in)
options = optimoptions('quadprog','Display','off');
persistent uk; %essa variavel pode reter valores entre "clear" do matlab. usar "clear all"
if(isempty(uk))
   uk = zeros(n_in,1); 
end

%% DEFINA AS RESTRICOES AQUI (2 entradas e 2 saidas)
Umax = [0.8;0.8]; %TITO
Umin = [-1;-1];
dUmax = [0.1;0.1];
dUmin = [-0.1;-0.1];
Ymax = [1; 1];
Ymin = [-0.5 ;-0.5];

% Umax = [5]; %SISO
% Umin = [0];
% Ymax = [1];
% Ymin = [-0.5];

M_inv = eye(nu*n_in);
n_ones = -1*ones(n_in*(nu-1),1);
n_diag = diag(n_ones,-n_in);
M_inv = M_inv+n_diag;
M = inv(M_inv);

u0 = [uk; zeros(n_in*(nu-1),1)];
Umax = repmat(Umax,nu,1);
Umin = repmat(Umin,nu,1);
Ymax = repmat(Ymax,n,1);
Ymin = repmat(Ymin,n,1);


%     %UMAX    
    A_con1 = [eye(nu*n_in);-eye(nu*n_in)];
    B_con1 = [Umax;-Umin];
%     %YMAX
%     A_con2 = [G;-G];
%     B_con2 = [Ymax-f;f-Ymin];
%     
%     %Aqui escolhemos quais das restrições anteriores são ustadas. Basta
%     %concatenar;
    A_con = [A_con1;[];[]];
    B_con = [B_con1;[];[]];
    
    U = quadprog(G'*Qd*G+M_inv'*Ql*M_inv,G'*Qd*(f-w)-M_inv'*Ql*u0,[],[],[],[],[],[],[],options); %sem restricoes
%      U = quadprog(G'*Qd*G+M_inv'*Ql*M_inv,G'*Qd*(f-w)-M_inv'*Ql*u0,A_con,B_con,[],[],[],[],[],options); %com restricoes

      % SATURADOR ATUA AQUI!
      if(U(1:n_in)> Umax(1:n_in))
          U(1:n_in) = Umax(1:n_in);
      elseif(U(1:n_in) < Umin(1:n_in))
          U(1:n_in) = Umin(1:n_in);
      end
      % FIM DO BLOCO SATURADOR
    
    uk = U(1:n_in);
end