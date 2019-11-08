clear all; close all; clc;

Ts=0.1;

global xo N Nu alfax alfay alfatheta Q R umin umax M ordemFilter
global ubase du ro xr yr tr vmax vmin wmax wmin delay_v delay_w alfaReal

% delays
delay_v = 0; delay_w = 0;

%%% Ajustes apenas para EPSAC. Ajustes João Emanuel estão
%%% na Sfunction MPC_SL, e ajustes do Kanklar estão na S-function
%%% Kanklar

% Horizontes de controle e predição
N = 5;
Nu = 1;

%%% Outros Ajustes
Qr = 0.01;
Qx = 1;
Qy = 1;
Qt = 0.01;

ordemFilter = 2;

%Initial Condition
xi = 0;
yi = 0;
thetai = 0;
xo = [xi;yi;thetai];

% saturation
vmax = 0.8;
vmin = -vmax;
wmax = 0.8;
wmin = -wmax;

% load reference
load ('ref.mat');
xr = ref(1,:);
yr = ref(2,:);
tr = ref(3,:);

% prepare data

%plots
% load reference
load ('ref.mat');
xr = ref(1,:);
yr = ref(2,:);
tr = ref(3,:);


%Motors PID
Kp=5;
Ki=10.0;
Kd=0;

% Parameters of the robot physic structure
r = 0.5/(2*pi); %raio da roda
d = 0.4; % comprimento do eixo entre rodas;

% Parameters of the DC motors and Speed Controller
J = 0.01; %kg.m^2  moment of inertia of the rotor
b = 0.1; %N.m.s   motor viscous friction constant
Ke = 0.08; % V/rad/sec  electromotive force constant
Kt = 0.08; % N.m/Amp   motor torque constant
Res = 1; % Ohm electric resistance
L = 0.5; %H  electric inductance

% weighting matrices
R = Qr*eye(2*Nu);
%R(1,1) = 0; R(2,2) = 0;
Q = diag([Qx Qy Qt]);
Qcell = repmat({Q},1,N);
Q = blkdiag(Qcell{:});

% load reference
load ('ref.mat');
load ('ubase.mat');
vr = ubase(1,:);
wr = ubase(2,:);
w = length(vr);
tmax = Ts*w-Ts;
t = 0:Ts:tmax;
k = 1:w;
k = [t' k'];
Tsim = tmax;

% account for extended prediction
for l = 1:N
    ref = [ref ref(:,end)];
    ubase = [ubase [0 0]'];
end
xr = ref(1,:);
yr = ref(2,:);
tr = ref(3,:);

% other definition

umin = [];
umax = [];

for ba = 1:Nu;
    umin = [umin;[vmin;wmin]];
    umax = [umax;[vmax;wmax]];
end

% M matrix
l = ones(1,2*Nu);
M = diag(l);

for na = 1:(2*Nu-2);
    M(na+2,na)=-1;
end;

% other definitions
du = 0.00001;

% limite de atualização
ro = 0.001;

angulos = 0:10:60;
angulos = angulos*pi/180;
alfaReal = -1.05;

figure(5);
plot(real(alfaReal),imag(alfaReal),'x')

for l = 1:length(angulos)

    alfa(l,1) = alfaReal;
    alfa(l,2) = alfaReal*exp(j*angulos(l));
    alfa(l,3) = alfaReal*exp(-j*angulos(l));
    
end

alfa = exp(alfa*Ts);


% legendas
legends = cell(length(alfa),1); % legends matrix

figure(1)
plot(xr(1:length(t)),yr(1:length(t)),'r--','LineWidth',2); hold on;

figure(2)
plot(t,tr(1:length(t)),'LineWidth',2); hold on;

figure(4)

subplot(2,1,1);
plot(t,vr(1:length(t)),'LineWidth',2); hold on;
subplot(2,1,2);
plot(t,wr(1:length(t)),'LineWidth',2); hold on;


for l = 1:length(alfa)
    
    % filter
    legends{l,1} = strcat('\theta = ',num2str(angulos(l)*180/pi),'°');
    alfax = alfa(l,2);
    alfay = alfa(l,2);
    alfatheta = alfa(l,2);
    alfaReal = alfa(l,1);
    
    % run simulation
    xo = [xi;yi;thetai];
    sim('motorsPID3');
    
    % rename simout data
    t = simout.time;
    v = simout.signals.values(:,1);
    w = simout.signals.values(:,2);
    
    % pos
    x = pos.signals.values(:,1);
    y = pos.signals.values(:,2);
    theta = pos.signals.values(:,3);
    
    % vel
    vreal = vel.signals.values(:,1);
    wreal = vel.signals.values(:,2);
    
    % Put theta in range [-pi,pi]
    for m = 1:length(theta)
        
        tmed = theta(m);
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
        
        theta(m) = tmed;
        
    end
    
    % prepare for plots
    
    figure(1)
    plot(x,y,'LineWidth',2); hold on;
    
    figure(2)
    plot(t,theta,'LineWidth',2); hold on;
    
    figure(3)
    err_x = x-xr(1:length(t))'; err_y= y-yr(1:length(t))';
    subplot(3,1,1);
    plot(t,err_x,'LineWidth',2); hold on;
    subplot(3,1,2);
    plot(t,err_y,'LineWidth',2); hold on;
    subplot(3,1,3);
    err_theta = atan2(sin(theta-tr(1:length(t))'),cos(theta-tr(1:length(t))'));
    plot(t,err_theta,'LineWidth',2); hold on;
    
    figure(4)
    
    subplot(2,1,1);
    v = [0;v(1:end-1)];
    stairs(t,v,'LineWidth',2); hold on;
    
    subplot(2,1,2);
    w = [0;w(1:end-1)];
    stairs(t,w,'LineWidth',2); hold on;
    
    if l > 1
        figure(5);
        plot([alfa(l,2) conj(alfa(l,2))],'x'); hold on;
    end
    
    Qe(l) = sum(err_x.^2) + sum(err_y.^2);
    disp(['\theta = ' num2str(angulos(l)*180/pi) ', Erro quadratico = ' num2str(Qe(l))])
    
end

figure(5)
zgrid

legends = ['reference'; legends];

figure(1)
ylabel('Y[m]');
xlabel('X[m]');
legend(legends);
axis tight;
grid on;

figure(2)
ylabel('Theta[rad]');
xlabel('T[s]');
legend(legends);
axis tight;
grid on;
%
figure(3)
subplot(3,1,1);
ylabel('X Error');
xlabel('T[s]');
legend(legends(2:end));
axis tight;
grid on;

subplot(3,1,2);
ylabel('Y Error');
xlabel('T[s]');
axis tight;
grid on;

subplot(3,1,3);
ylabel('Theta Error');
xlabel('T[s]');
axis tight;
grid on;

%
figure(4)
subplot(2,1,1);
legend(legends);
axis tight;
grid on;
subplot(2,1,2);
xlabel('T[s]');
axis tight;
grid on;

figure(5)

Qe = sum(err_x.^2) + sum(err_y.^2) + sum(err_theta.^2);