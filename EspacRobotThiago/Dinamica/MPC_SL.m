function [sys, x0, str, ts] = MPC_SL(t,x,u,flag,ts)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%inicialização
if flag == 0
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 2;
    sizes.NumInputs      = 4;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    
    sys = simsizes(sizes);
    
    x0 = [ ];
    str = [ ];
    ts = [0.1]; %tempo de amostragem variável
    
    %Calcula próximo instante de amostragem
elseif flag == 4
    sys=[];
    
elseif flag == 3
    
    % other global variables
    global u0_MPC N ubase xr yr
    global tr vmax vmin wmax wmin
    
    % interaction control
    k = round(u(4));
    
    % Put theta in range [-pi,pi]
    tmed = u(3);
    Npi = tmed/(2*pi);
    
    if Npi > 1;
        tmed = tmed - (round(Npi))*2*pi;
    end
    
    if Npi < -1;
        tmed = tmed - (round(Npi))*2*pi;
    end
    
    if tmed > pi
        tmed = tmed - 2*pi;
    end
    
    if tmed < -pi
        tmed = tmed + 2*pi;
    end
    
    y = [u(1);u(2);tmed]; %% medição, entrada de parâmetro
    
    lb = [];
    ub = [];
    
    Qe = diag([1 1 0.1]);
    Re = 0.1*eye(2);
    
    for l = 1:N
     Qb(3*l-2:3*l,3*l-2:3*l) =Qe;
     Rb(2*l-1:2*l,2*l-1:2*l) =Re;
     lb = [lb;[vmin;wmin]];
     ub = [ub;[vmax;wmax]];
    end
    
    A = [];
    xb = [];
    aux = eye(3);
    
    Ts = 0.1;
    for j = 0:N-1
 

        Aj{j+1} = [1 0 -Ts*ubase(1,k+j)*sin(tr(k+j)); 0 1  Ts*ubase(1,k+j)*cos(tr(k+j)); 0 0  1];
        
        Bj{j+1} = [cos(tr(k+j))*Ts  0;
            sin(tr(k+j))*Ts  0;
            0                Ts];
        aux = Aj{j+1}*aux;
        A = [A;aux];
        xb = [xb; ubase(:,k+j)]; %************************
    end
    
    
    B = [];
    
    for j = 1:N
        aux2  = eye(3);
        Baux = [];
        for i = 1:j
            if i<j
                Baux = [aux2*Bj{j-i+1},Baux];
                aux2 = aux2*Aj{j-i+1};
            else
                Baux = [aux2*Bj{j-i+1},Baux,zeros(3,(N-j)*2)];
            end
        end
        B = [B;Baux];
    end
    
    H = 2*(B'*Qb*B+Rb);
    
    Xr = [xr(k);yr(k);tr(k)];
    
    % reference errors
    errx = y(1)-xr(k);
    erry = y(2)-yr(k);
    errt = y(3)-tr(k);
    errt = atan2(sin(errt), cos(errt)); % smooth the error
    E = [errx;erry;errt];
    
    f = 2*B'*Qb*A*E;
    options = optimset('Display','off','Largescale','off');
    Uopt = quadprog(H,f,[],[],[],[],[lb-xb],[ub-xb],[],options);
    %Uopt = -inv(H)*f;
    
    u = ubase(:,k)+Uopt(1:2,1);
    
    % saturation
    u(1,1) = min(u(1,1),vmax);
    u(1,1) = max(u(1,1),vmin);
    u(2,1) = min(u(2,1),wmax);
    u(2,1) = max(u(2,1),wmin);
    
    % control action
    u0_MPC = u(:,1);
    sys = u0_MPC;
    
else
    sys = [ ]; %não faz nada
end

end

