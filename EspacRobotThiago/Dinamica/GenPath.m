clear all; close all; clc;

global Ts;
Ts = 0.1;
xo = [0 0 0]'; % posi��o inicial
ref(:,1) = [0 0 0]';
disp('Escolha a trajet�ria de acordo com as op��es abaixo: '); 
traj = input('  1 - circular \n  2 - oito\n  3 - em L\n  4 - Quadrado\n  5 - Sen�ide\n\n Digite o valor correspondente � trajet�ria desejada: ');

N = 5;
% velocity profile
v = 0.15;
w = 0.15;
u1 = [v w]'; % for circle
u2 = [v -w]'; % for circle and eight
u3 = [v 0]'; % for L and square

m = input('Digite o n�mero de voltas/per�odos desejadas: '); %number of laps

if (traj == 2) || (traj == 1);  % eight or circle
    k = round(traj*(2*pi/(abs(u1(2))*Ts)+1));
    
elseif traj == 3
    % L trajectory
    L1 = input('Entre o comprimento do primeiro trecho: '); %[m]
    L2 = input('Entre o comprimento do segundo trecho: '); %[m]
    k1 = (L1/v)/Ts;
    k2 = (L2/v)/Ts;
    k = k1+k2;
    
elseif traj == 4 % square
    L = input('Entre o comprimento do lado: '); %[m]
    k1 = (L/v)/Ts; k2 = 2*k1; k3 = 3*k1; k4 = 4*k1;
    k = 4*k1;

elseif traj == 5 %senoidal
    A = input('Entre a amplitude da senoide: ');
    L = input('Entre o comprimento m�ximo: ');
    b = m*(2.0*pi/L);
    k = round(L/(w*Ts))+1;
end

t = 0:Ts:Ts*k-Ts+5*N*Ts;

for l = 2:(k+5*N); 
    
    % para oito
    if traj == 2;
        
        if l >= k/2
            p = modelo(u2,xo,Ts);
            ubase(:,l) = u2;
        else
            p = modelo(u1,xo,Ts);
            ubase(:,l) = u1;
        end
      
    % para circulo
    elseif traj == 1;
        p = modelo(u1,xo,Ts);
        ubase(:,l) = u1;
    
    % para L
    elseif traj == 3;
        if l<=k1;
            xo(3) = 0;
            p = modelo(u3,xo,Ts);
            ubase(:,l) = u3;
        else
            xo(3) = pi/2;
            p = modelo(u3,xo,Ts);
            ubase(:,l) = u3;
        end
    
    % para quadrado
    elseif traj == 4;
        if l <=k1;
            xo(3) = 0;
            p = modelo(u3,xo,Ts);
            ubase(:,l) = u3;
        elseif l >k1 && l <=k2;
            xo(3) = pi/2;
            p = modelo(u3,xo,Ts);
            ubase(:,l) = u3;
        elseif l >k2 && l <=k3;
            xo(3) = -pi;
            p = modelo(u3,xo,Ts);
            ubase(:,l) = u3;
        elseif l >k3;
            xo(3) = -pi/2;
            p = modelo(u3,xo);
            ubase(:,l) = u3;
        end
        
    % para sen�ide
    elseif traj == 5;
   
        x(l) = w*t(l);
        y(l) = A*sin(b*w*t(l));
        Xr_dot(l) = w;
        Yr_dot(l) = A*b*w*cos(b*w*t(l));
        Xr_dot2(l) = 0;
        Yr_dot2(l) = -A*b^2*w^2*sin(b*w*t(l));
        theta(l) = atan2(Yr_dot(l), Xr_dot(l));
        Vr(l) = sqrt(Xr_dot(l)^2 + Yr_dot(l)^2);
        Wr(l) = (Yr_dot2(l)*Xr_dot(l)-Xr_dot2(l)*Yr_dot(l))/(Xr_dot(l)^2+Yr_dot(l)^2);
        ubase(:,l) = [Vr(l);Wr(l)];
        p = [x(l);y(l);theta(l)];
    end
    

    x(l)=p(1);
    y(l) = p(2);
    theta(l) = p(3);
    xo = p;
    % save output references
    ref(:,l) = p;

end;

ref = repmat(ref,1,m); % account for the number of laps
ubase = repmat(ubase,1,m);

save('ref.mat','ref');
save('ubase.mat','ubase');

figure(1)
subplot (1,3,1);
stairs(x,y,'r','LineWidth', 2);
legend('Reference Trajectory');
ylabel('Y[m]');
xlabel('X[m]');
axis tight;
grid on;

subplot (1,3,2);
plot(t,theta(1:length(t)),'LineWidth', 2);
legend('Theta');
ylabel('Theta [rad]');
axis tight;
grid on;

subplot (1,3,3);
plot(t,ubase(1,1:length(t)),t,ubase(2,1:length(t)),'LineWidth', 2);
legend('Vr','Wr');
ylabel('Vr[m/s],Wr[rad/s]');
axis tight;
grid on;