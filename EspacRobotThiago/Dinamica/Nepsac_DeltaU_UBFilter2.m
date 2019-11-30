function [sys, x0, str, ts] = Nepsac_DeltaU_UBFilter2(t,x,u,flag,ts)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%inicializa��o
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
    ts = [0.1]; %tempo de amostragem vari�vel
    
    %Calcula pr�ximo instante de amostragem
elseif flag == 4
    sys=[];
    
elseif flag == 3
    
    % filter global variables
    global zix ziy zitheta Cx Cy Ctheta npx npy nptheta
    global Fx Fy Ftheta alfax alfay alfatheta alfaReal
    
    % other global variables
    global xo u0_EPSAC N Nu Q R M ubase du ro xr yr Ts
    global tr vmax vmin wmax wmin delay_v delay_w ordemFilter
    
    % interaction control
    k = round(u(4));
    
    % initializes filter and other parameters
    if k == 1;
        u0_EPSAC = [0;0];
        % disturbance filtering x
        D = [1 -1 zeros(1,ordemFilter-1)];
        
        if (ordemFilter == 1);
            
            Cx = [1 -alfax];
            
        elseif (ordemFilter == 2);
            Cx = conv([1 -alfax],[1 -conj(alfax)]);
            
        elseif (ordemFilter == 3);
            Cx = conv([1 -alfaReal],conv([1 -alfax],[1 -conj(alfax)]));
            
        end
        
        [~,Fx] = GetFE(Cx,D,5*N);
        Fx = [Fx(1,:); Fx(3,:); Fx(6,:); Fx(10,:); Fx(15,:)]
        npx = zeros(length(Cx)-1,1);
        [~,zix] = filter(1,Cx,0);
        
        % disturbance filtering y
        
        if (ordemFilter == 1);
            Cy = [1 -alfay];
        elseif (ordemFilter == 2);
            Cy = conv([1 -alfay],[1 -conj(alfay)]);
        elseif (ordemFilter == 3);
            Cy = conv([1 -alfaReal],conv([1 -alfay],[1 -conj(alfay)]));
        end
        
        [~,Fy] = GetFE(Cy,D,5*N);
        Fy = [Fy(1,:); Fy(3,:); Fy(6,:); Fy(10,:); Fy(15,:)]
        npy = zeros(length(Cy)-1,1);
        [~,ziy] = filter(1,Cy,0);
        
        % disturbance filtering theta
        if (ordemFilter == 1);
            Ctheta = [1 -alfatheta];
        elseif (ordemFilter == 2);
            Ctheta = conv([1 -alfatheta],[1 -conj(alfatheta)]);
        elseif (ordemFilter == 3);
            Ctheta = conv([1 -alfaReal],conv([1 -alfatheta],[1 -conj(alfatheta)]));
        end
        
        [~,Ftheta] = GetFE(Ctheta,D,5*N);
        Ftheta = [Ftheta(1,:); Ftheta(3,:); Ftheta(6,:); Ftheta(10,:); Ftheta(15,:)]
        nptheta = zeros(length(Ctheta)-1,1);
        [~,zitheta] = filter(1,Ctheta,0);
        
    end
    
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
    
    y = [u(1);u(2);tmed]; %% medi��o, entrada de par�metro
    
    %nada haver com o controlador. Pr�ximo passo � determinar resposta base

    
    x=modelo(u0_EPSAC,xo,0.1); % (CALCULADA)
    n = y - x;
    xo = y;
    
    % filter x
    [nfx,zfx] = filter(1,Cx,n(1),zix);
    zix = zfx;
    npx = [nfx;npx(1:end-1)];
    Nx = Fx*npx;
    
    % filter y
    [nfy,zfy] = filter(1,Cy,n(2),ziy);
    ziy = zfy;
    npy = [nfy;npy(1:end-1)];
    Ny = Fy*npy;
    
    % filter theta
    [nftheta,zftheta] = filter(1,Ctheta,n(3),zitheta);
    zitheta = zftheta;
    nptheta = [nftheta;nptheta(1:end-1)];
    Ntheta = Ftheta*nptheta;
    
    nfiltro = [Nx'; Ny'; Ntheta'];
    
    % Base update
    if k > 1
        urefp = ubase(:,k-1);
    end
    
    uref = ubase(:,k:N+k-1);
    
    ub = zeros(2,N);
    
    ub(1,:) = uref(1,1); % lembrar de verificar
    ub(2,:) = 0;
    
    vbd = ub;
    vbd(1,1) = vbd(1,1) + du;
    
    wbd = ub;
    wbd(2,1) = wbd(2,1) + du;
    
    % filtros
    
    ubf = delayFilter(ub,delay_v,delay_w);
    vbdf = delayFilter(vbd,delay_v,delay_w);
    wbdf = delayFilter(wbd,delay_v,delay_w);
    %     ubf = ubfilter(ub);
    %     vbdf = ubfilter(vbd);
    %     wbdf = ubfilter(wbd);
    
    Yb = ybase(ubf,xo,nfiltro);
    xb = Yb(1,:);
    yb = Yb(2,:);
    tb = Yb(3,:);
    
    % counter
    aux = ro+0.1*ro;
    
    while (aux > ro)
        
        % G matrix
        % impulse response
        xbv = ybase(ubf,xo,nfiltro); xbv = reshape(xbv,[3*N,1]);
        xbvd = ybase(vbdf,xo,nfiltro); xbvd = reshape(xbvd,[3*N,1]);
        xbw = ybase(ubf,xo,nfiltro); xbw = reshape(xbw,[3*N,1]);
        xbwd = ybase(wbdf,xo,nfiltro); xbwd = reshape(xbwd,[3*N,1]);
        
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
        xerro = [xr(k+1) xr(k+3) xr(k+6) xr(k+10) xr(k+15)];
        yerro = [yr(k+1) yr(k+3) yr(k+6) yr(k+10) yr(k+15)];
        thetaerro = [tr(k+1) tr(k+3) tr(k+6) tr(k+10) tr(k+15)];
%         errx = (xb-xr(:,k:N+k-1))'; % n deveria ser k+1:N+k? checar
%         erry = (yb-yr(:,k:N+k-1))';
%         errt = (tb-tr(:,k:N+k-1))';
% xb
% xerro
        errx = (xb-xerro)'; % n deveria ser k+1:N+k? checar
        erry = (yb-yerro)';
        errt = (tb-thetaerro)';

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
            L = [urefp-u0_EPSAC; zeros(2*Nu,1)];
        else
            L = [-u0_EPSAC; zeros(2*Nu,1)];
        end
        
        L = L(1:2*Nu,:); % account for difference btw Nu and N
        Ku = M*UbUr + L;    % mgn
        Fnep = 2*(Gnep'*Q*E + M'*R*Ku);   % mgn
        
%         % trial
% %         original cost function
%                 Hnep = 2*(Gnep'*Q*Gnep+R);
%                 Fnep = 2*(Gnep'*Q*E+R*UbUr);
 
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
        aux = sum(uo*uo');
        aux = 0; % = 0 -> epsac
    end
    
    u = ub(:,1);
    
    % saturation
    u(1,1) = min(u(1,1),vmax);
    u(1,1) = max(u(1,1),vmin);
    u(2,1) = min(u(2,1),wmax);
    u(2,1) = max(u(2,1),wmin);
    
    % control action
    u0_EPSAC = u(:,1);
    sys = u0_EPSAC;
    
else
    sys = [ ]; %n�o faz nada
end

end

