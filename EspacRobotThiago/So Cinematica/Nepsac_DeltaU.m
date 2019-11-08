clear all; close all; clc;
global Ts;

N = 1; % Horizonte de predi��o
Nu = 1; % Horizonte de controle

% initial pose
xo = [0; -0.5; pi/2];

% saturation
vmax = 0.8;
vmin = 0;
wmax = 0.8;
wmin = -wmax;

umin = [];
umax = [];

for k = 1:Nu;
    umin = [umin;[vmin;wmin]];
    umax = [umax;[vmax;wmax]];
end

% weighting matrices
R = 0.01*eye(2*Nu);
%R(1,1) = 0; R(2,2) = 0;
Q = diag([1 1 0.01]);
Qcell = repmat({Q},1,N);
Q = blkdiag(Qcell{:});

% M matrix
l = ones(1,2*Nu);
M = diag(l);

for na = 1:(2*Nu-2);
    M(na+2,na)=-1;
end;

% load reference
load ('ref.mat');
w = length(ref);
load ('ubase.mat');

% account for extended prediction
for l = 1:N
    ref = [ref ref(:,end)];
    ubase = [ubase [0 0]'];
end
xr = ref(1,:); yr = ref(2,:); tr = ref(3,:);

% other definitions
du = 0.00001;
ub = ubase(:,1:N);% tamanho do horizonte de predi��o
ub(1,:) = 0.15;
ub(2,:) = 0;
%u0 = ub(:,1);
u0 = [0 0]';

% time profile
Ts = 0.1; 
t = 0:Ts:Ts*w-Ts;

% control disturbs
Tpert = 50;
Psamples = 1; % # samples with disturb
P = [0 0]';

% limite de atualiza��o
ro = 0.001;

% Creates noise profile
Mean = 0; % zero mean
sd_xy = 0.001; % standard deviation
sd_t = 0.001; % standard deviation
noise_xy = Mean + sd_xy.*randn(2,w);
noise_t = Mean + sd_t.*randn(1,w);
noise = [noise_xy;noise_t];

noise = 0*noise; % no noise

for k = 1:w;
    
    if k >= Tpert/Ts && k <= (Tpert/Ts + Psamples)
        y(:,k)=modelo(u0+P,xo)+noise(:,k); % *(LIDA) simulando a sa�da da planta,
    else
        y(:,k)=modelo(u0,xo)+noise(:,k);
    end

    %nada haver com o controlador. Pr�ximo passo � determinar resposta base 
    x(:,k)=modelo(u0,xo); % (CALCULADA)
    n = y(:,k) - x(:,k);
    
%     [nf(k),zf] = filter(1,C,n(k),zi);
%     zi = zf;
%     np = [nf(k);np(1:end-1)];
%     Np = F*np;
    
    xo = y(:,k);
    
    % Base update
    if k > 1
        urefp = uref(:,1);
    end
    
    uref = ubase(:,k:N+k-1);
    ub = [ub(:,2:end) ub(:,end)]; %ignorei no python
    
    ub(1,:) = 0.3; % lembrar de verificar
    ub(2,:) = 0; % lembrar de verificar
    
    vbd = ub;
    vbd(1,1) = vbd(1,1) + du;
    wbd = ub;
    wbd(2,1) = wbd(2,1) + du;
    
    Yb = ybase(ub,xo,n); 
    xb = Yb(1,:);
    yb = Yb(2,:);
    tb = Yb(3,:);
    
    % counter
    co = 0;
    aux = ro+0.1*ro;
    
    while (aux > ro)
            
        % G matrix
        % impulse response
        xbv = ybase(ub,xo,n); xbv = reshape(xbv,[3*N,1]);
        xbvd = ybase(vbd,xo,n); xbvd = reshape(xbvd,[3*N,1]);
        xbw = ybase(ub,xo,n); xbw = reshape(xbw,[3*N,1]);
        xbwd = ybase(wbd,xo,n); xbwd = reshape(xbwd,[3*N,1]);
   
        % for v
        gv = (xbvd-xbv)/du; 
        
        % for w
        gw = (xbwd-xbw)/du; 
            
        %
        Gnep = zeros(3*N,2*N);
        nep = [];
        av = reshape(gv,[3,N]);
        aw = reshape(gw,[3,N]);
        
        for j = 1:N
            nep = [nep;[av(:,j),aw(:,j)]];
            Gnep((3*(N-j)+1):end,2*(N-j)+1:2*(N-j+1)) = nep;
        end
        
        % for Nu different than N
        if (N ~= Nu);
            % step response
            avg = av(:,1); awg = aw(:,1);
            for l = 2:N
                avg = [avg avg(:,l-1)+av(:,l)];
                awg = [awg awg(:,l-1)+aw(:,l)];
            end
            avg = reshape(avg,[3*N,1]);
            awg = reshape(awg,[3*N,1]);
            Gnep = Gnep(:,1:2*Nu);
            Gnep(:,end-1) = [zeros(3*(Nu-1),1);avg(1:end-3*(Nu-1))];
            Gnep(:,end) = [zeros(3*(Nu-1),1);awg(1:end-3*(Nu-1))];
        end
        
        % reference errors
        errx = (xb-xr(:,k:N+k-1))';
        erry = (yb-yr(:,k:N+k-1))';
        errt = (tb-tr(:,k:N+k-1))';
        errt = atan2(sin(errt), cos(errt)); % smooth the error
    
        E = [];
        for s = 1:length(errx)
            E = [E errx(s) erry(s) errt(s)];
        end
        E = E';

        % cost function
        
        Hnep = 2*(Gnep'*Q*Gnep+M'*R*M);
        xbnep = reshape(ub(:,1:Nu),[2*Nu,1]);
        UbUr = xbnep-reshape(uref(:,1:Nu),[2*Nu,1]);
        
        if k>1
        L = [urefp-u0; zeros(2*Nu,1)];               
        else
        L = [-u0; zeros(2*Nu,1)];                
        end
        
        L = L(1:2*Nu,:); % account for difference btw Nu and N
        Ku = M*UbUr + L;    % mgn
        Fnep = 2*(Gnep'*Q*E + M'*R*Ku);   % mgn
        
        % analytical solution
        uo = -Hnep\Fnep;
        
        % quadprog solution
        %options = optimset('Display','off','Largescale','off');
        %uo = quadprog(Hnep,Fnep,[],[],[],[],[umin-xbnep],[umax-xbnep],[],options);
    
        % account for Nu < N
        uo = reshape(uo,[2,Nu]);
        [muo,nuo] = size(uo); [mub,nub] = size(ub); na = nub - nuo;
        uo = [uo [uo(1,end).*ones(1,na);uo(2,end).*ones(1,na)]]; %for when N ~= Nu
        
        % update Ubase and counter
        ub = uo+ub;
        co = co + 1;
        aux = sum(uo*uo');
        aux = 0; % epsac
    end
    
    iteracoes(k) = co;
    u(:,k) = ub(:,1);
    
    % saturation
    u(1,k) = min(u(1,k),vmax);
    u(1,k) = max(u(1,k),vmin);
    u(2,k) = min(u(2,k),wmax);
    u(2,k) = max(u(2,k),wmin);
    
    % control action
    u0 = u(:,k); 

end

[~,n] = size(y);
%MSE(mm)=sqrt((sum(sum((y-ref(:,1:n)).^2)))/(3*n)); 
MSE=sum(sum((y(1:2,:)-ref(1:2,1:n)).^2));

%plots
figure(1)
%subplot(2,2,1)
plot(y(1,:),y(2,:),'k','LineWidth',2); hold on;
plot(ref(1,1:length(y(1,:))),ref(2,1:length(y(1,:))),'r--','LineWidth',2);
legend('Real Robot','Reference Robot');
ylabel('Y[m]');
xlabel('X[m]');
axis tight;
grid on;

% subplot(2,2,3)
% plot(t(1:k),y(3,:),'LineWidth',2); hold on;
% plot(t(1:k),tr(1:k),'LineWidth',2); 
% legend('Real Robot','Reference Robot');
% ylabel('Theta[rad]');
% xlabel('T[s]');
% axis tight;
% grid on;
% 
figure(2)
plot(t(1:k),y(1,:)-xr(1:length(y)),'LineWidth',2); hold on;
plot(t(1:k),y(2,:)-yr(1:length(y)),'LineWidth',2); hold on;
err_theta = atan2(sin(y(3,:)-tr(1:length(y))),cos(y(3,:)-tr(1:length(y))));
plot(t(1:k),err_theta,'LineWidth',2); 
legend('Xerro','Yerro','Thetaerro');
ylabel('Error');
xlabel('T[s]');
axis tight;
grid on;
% 
% subplot(2,2,4)
figure(3)
stairs(t(1:k),ubase(1,1:k),'LineWidth',2); hold on;
stairs(t(1:k),u(1,1:k),'LineWidth',2); hold on;
stairs(t(1:k),ubase(2,1:k),'LineWidth',2); hold on;
stairs(t(1:k),u(2,1:k),'LineWidth',2);
legend('vref','v','wref','w');
xlabel('T[s]');
axis tight;
grid on;

% figure(3);
% stairs(t,iteracoes,'LineWidth', 2);
% legend('Interactions');
% ylabel('# interactions');
% axis tight;
% grid on;

disp(['O erro quadr�tico para cada trajet�ria �:']);
disp(num2str(MSE));
disp(['O erro quadr�tico para todas as trajet�rias �:']);
disp(num2str(sum(MSE)));