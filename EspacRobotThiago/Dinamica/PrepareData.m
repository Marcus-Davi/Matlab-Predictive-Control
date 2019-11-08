
% EPSAC
% control
t = con_EPSAC.time;
v_EPSAC = con_EPSAC.signals.values(:,1);
w_EPSAC = con_EPSAC.signals.values(:,2);

% pose
x_EPSAC = pos_EPSAC.signals.values(:,1);
err_x_EPSAC = x_EPSAC-xr(1:end-5*N)'; 
y_EPSAC = pos_EPSAC.signals.values(:,2);
err_y_EPSAC= y_EPSAC-yr(1:end-5*N)';
theta_EPSAC = pos_EPSAC.signals.values(:,3);

% Put theta in range [-pi,pi]
for l = 1:(length(theta_EPSAC))
    
    tmed = theta_EPSAC(l);
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
    
    theta_EPSAC(l) = tmed;
    
end

err_theta_EPSAC = atan2(sin(theta_EPSAC-tr(1:end-5*N)'),cos(theta_EPSAC-tr(1:end-5*N)'));

% Quadratic error info

Qe_EPSAC = costF(err_x_EPSAC,err_y_EPSAC,err_theta_EPSAC);

disp(['Quadratic Error for EPSAC was ',num2str(Qe_EPSAC)]);

plota