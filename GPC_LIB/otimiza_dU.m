function dU = otimiza_dU(G,Ql,Qd,f,w,nu,n,n_in)
options = optimoptions('quadprog','Display','off');
persistent uk; %essa variavel pode reter valores entre "clear" do matlab. usar "clear all"
if(isempty(uk))
   uk = zeros(n_in,1); 
end

%% DEFINA AS RESTRICOES AQUI (2 entradas e 2 saidas)
Umax = [0.8;0.2]; %TITO
Umin = [-1;-1];
dUmax = [0.1;0.1];
dUmin = [-0.1;-0.1];
Ymax = [1; 1];
Ymin = [-0.5 ;-0.5];
% 
% Umax = [0.4]; %SISO
% Umin = [0];
% dUmax = [0.1];
% dUmin = [-0.1];
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
dUmax = repmat(dUmax,nu,1);
dUmin = repmat(dUmin,nu,1);
Ymax = repmat(Ymax,n,1);
Ymin = repmat(Ymin,n,1);

    A_con1 = [M;-M];
    B_con1 = [Umax-u0;u0-Umin];
    %DUMAX    
    A_con2 = [eye(nu*n_in);-eye(nu*n_in)];
    B_con2 = [dUmax;-dUmin];
    %YMAX
    A_con3 = [G;-G];
    B_con3 = [Ymax-f;f-Ymin];
    
    %Aqui escolhemos quais das restrições anteriores são ustadas. Basta
    %concatenar;
    A_con = [A_con1;[];[]];
    B_con = [B_con1;[];[]];
    
%     dU = quadprog(G'*G+Ql,G'*(f-w),[],[],[],[],[],[],[],options); %sem restricoes
    dU = quadprog(G'*Qd*G+Ql,G'*Qd'*(f-w),[],[],[],[],[],[],[],options); %com restricoes
    
    uk = uk + dU(1:n_in);
    

end