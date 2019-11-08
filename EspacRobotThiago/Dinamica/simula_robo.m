clear all; close all; clc;

Ts=0.1;

global xo N Nu alfax alfay alfatheta Q R umin umax M ordemFilter
global ubase du ro xr yr tr vmax vmin wmax wmin delay_v delay_w

% delays
delay_v = 0; delay_w = 0;

%%% Ajustes apenas para EPSAC. Ajustes Jo�o Emanuel est�o
%%% na Sfunction MPC_SL, e ajustes do Kanklar est�o na S-function
%%% Kanklar

% Horizontes de controle e predi��o
N = 5;
Nu = 1;

%%% Outros Ajustes
Qr = 0.001;
Qx = 10;
Qy = 10;
Qt = 0.00;

ordemFilter = 2;

%Initial Condition
xi = 0;
yi = -0.5;
thetai = pi/2;
xo = [xi;yi;thetai];

% saturation
vmax = 0.8;
vmin = -vmax;
wmax = 0.8;
wmin = -wmax;

% filter
alfax = 0.9+0*1i;
alfay = 0.9+0*1i;
alfatheta = 0.9+0*1i;

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
w = length(vr)-5*N;
tmax = Ts*w-Ts;
t = 0:Ts:tmax;
k = 1:w;
k = [t' k'];
Tsim = tmax;

% % account for extended prediction
% for l = 1:N
%     %ref = [ref ref(:,end)];
%     %ubase = [ubase [0 0]'];
%     ubase = [ubase ubase(:,end)];
% end

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

% limite de atualiza��o
ro = 0.001;

% run simulation
sim('motorsPID2');

% prepare for plots

% load reference
load ('ref.mat');
xr = ref(1,:);
yr = ref(2,:);
tr = ref(3,:);

% prepare data

PrepareData