function [sys, x0, str, ts] = Kanklar(t,x,u,flag,ts)
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
    global u0 ubase xr yr
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
    
    % medição, entrada de parâmetro
    
    x = u(1); y = u(2); theta = tmed;
    
    % Kanklar algorithm
    
    % calcula erros
    
    e_x = xr(:,k) - x;
    e_y = yr(:,k) - y;
    e_theta = tr(:,k) - theta;
    e_theta = atan2(sin(e_theta), cos(e_theta)); % smooth the error
    
    % mudança de frame
    
    e1 = cos(theta)*e_x+sin(theta)*e_y;
    e2 = cos(theta)*e_y-sin(theta)*e_x;
    e3 = e_theta;
    
    % lei de controle
    
    % parâmetros de ajuste
  
    ksi = 0.8; % ajuste otimo 0.8
    g = 30; % ajuste ótimo 30
    
    % ganhos
    ref_v = ubase(1,k);
    ref_w = ubase(2,k);
    
    w_n = sqrt(ref_w^2 + g*ref_v^2);
    k1 = 2*ksi*w_n;
    k3 = k1;
    
    k2 = g*abs(ref_v);
    
    v1 = -k1*e1;
    v2 = -sign(ref_v)*k2*e2-k3*e3;
    
    vu = ref_v*cos(e3)-v1;
    wu = ref_w-v2;
    
    u = [vu;wu];
    
    % saturation
    u(1,1) = min(u(1,1),vmax);
    u(1,1) = max(u(1,1),vmin);
    u(2,1) = min(u(2,1),wmax);
    u(2,1) = max(u(2,1),wmin);
    
    % control action
    u0 = u(:,1);
    sys = u0;
    
else
    sys = [ ]; %não faz nada
end

end

