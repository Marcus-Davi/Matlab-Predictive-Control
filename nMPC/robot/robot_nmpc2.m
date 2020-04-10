clear;close all;clc
%% SIMULATION PARAMETERS
Ts = 0.1;
v = 0.1;
% R = 1;
% [Xr,Ur,Tsim] = path_oito(R,v,Ts);
[Xr,Ur,Tsim] = path_U(6.5,v,Ts,[0 0 pi/2]');
iterations = round((Tsim/Ts));

% plot(Xr(1,:),Xr(2,:),'k','linewidth',2)
% return
%% Controller Parameters
N = 5;
Nu = 5;
n_in = 2; %[v,w]
n_out = 3; %[x,y,theta]

Qx = 1;
Qy = 1;
Qt =  0.05;
Q = diag([Qx Qy Qt]);
for i=1:N
Qcell{i} = 2^(1-1)*Q;
Qncell{i} = 2^(N-1)*Q;
end
Ql = blkdiag(Qcell{:}); %Ref
% return
Qv = 0.1;
Qw = 0.1;
P = 0*blkdiag(Qncell{:});


Qu = diag([Qv Qw]);
Qcell = repmat({Qu},1,Nu);
Qu = blkdiag(Qcell{:}); %Ref

% Ql
% return

%% Simulation
yk = [-1 -1 0]';
u = [0. 0]';
YK = zeros(n_out,iterations);
umax = [0.47 3.77]';
opts = optimoptions('fmincon','Display','off');

u0 = repmat(u,N,1);
%Ax < b

D0 = [eye(n_in);-eye(n_in)]
D = []
for i=1:N
   D = blkdiag(D,D0) ;
end

b = [umax;umax];
b = repmat(b,N,1);

for k=1:iterations
    
       
    %Future References
    [Wr,Ur_] = getRef(Xr,Ur,k,N,Nu);
    C =  @(u) cost(Wr,yk,Ur_,u,N,Nu,Ql,Qu,P,Ts);
    
    J = cost(Wr,yk,Ur_,u0,N,Nu,Ql,Qu,P,Ts)
    
    %  Gauss - Newton

%     gn_i = 20;
%     u0 = repmat(u,N,1);
%     
%     for j = 1:gn_i
%         
%         Yp = [];
%         yp = yk;
%         u0_cell = mat2cell(u0,n_in*ones(1,Nu));
% %         celldisp(u0_cell);
%         for i=1:N
%            yp = robot_model(yp,u0_cell{i},Ts); %prediction k+i
%            Yp = [Yp;yp];
%            jac{i}  = [Ts*cos(yp(3)) 0;Ts*sin(yp(3)) 0;0 Ts];
%            
%         end
%         Jac = Ql*blkdiag(jac{:});
%         
%         u0 = u0 - inv(Jac'*Jac)*Jac'*Ql*(Yp - Wr);
%         J = (Yp-Wr)'*Ql*(Yp-Wr)
%     end
           
    
     % fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
   u0 = fmincon(C,u0,D,b,[],[],[],[],[],opts);

    u = u0(1:n_in);
    
%     if u(1) > 0.5
%         u(1) = 0.5;
%     end
%     
%     if(u(2) > 1)
%         u(2) = 1;
%     elseif u(2) < -1
%         u(2) = -1;
%     end
    yk = robot_model(yk,u,Ts); %yk = measurement
    
    
    YK(:,k) = yk;
    UK(:,k) = u;
    
end
%% Plots


plot(YK(1,1:k),YK(2,1:k),'b','linewidth',1)
hold on
plot(Xr(1,1:k),Xr(2,1:k),'k--','linewidth',2)
figure
grid on
plot(UK(1,:))
hold on
plot(UK(2,:))

%% Utils

function [wr,ur] = getRef(Xr,Ur,k,N,Nu)
wr = [];
ur = [];
[~,kend] = size(Xr);
for i=1:N    
     %testa final da traj. se sim, repete
    if(k+i < kend)
    wr = [wr;Xr(:,k+i)];    
%     ur = [ur;Ur(:,k+i)];    
    else
    wr = [wr;Xr(:,end)];
%     ur = [ur;Ur(:,end)];
    end

end

for i=1:Nu
     %testa final da traj. se sim, repete
    if(k+i < kend)
%     wr = [wr;Xr(:,k+i)];    
    ur = [ur;Ur(:,k+i)];    
    else
%     wr = [wr;Xr(:,end)];
    ur = [ur;Ur(:,end)];
    end

end

end

function y = cost(Wr,yk,Ur,u,N,Nu,Ql,Qu,P,Ts)
    %Predictions
    Yp = [];
    yp = yk;
    
    for i=1:N
       yp = robot_model(yp,u(2*(i-1)+1:2*(i-1)+2),Ts); %prediction k+i
       Yp = [Yp;yp];
    end
    
    

y = (Yp-Wr)'*Ql*(Yp-Wr) + (u-Ur)'*Qu*(u-Ur) + (Yp-Wr)'*P*(Yp-Wr);


end
