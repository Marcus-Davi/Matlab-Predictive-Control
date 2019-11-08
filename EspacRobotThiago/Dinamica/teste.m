% MPC
% control
v_MPC = con_MPC.signals.values(:,1);
w_MPC = con_MPC.signals.values(:,2);

% pose
x_MPC = pos_MPC.signals.values(:,1);
err_x_MPC = x_MPC-xr'; 
y_MPC = pos_MPC.signals.values(:,2);
err_y_MPC= y_MPC-yr';
theta_MPC = pos_MPC.signals.values(:,3);

% Put theta in range [-pi,pi]
for l = 1:length(theta_MPC)
    
    tmed = theta_MPC(l);
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
    
    theta_MPC(l) = tmed;
    
end

err_theta_MPC = atan2(sin(theta_MPC-tr'),cos(theta_MPC-tr'));

% Quadratic error info

Qe_MPC = costF(err_x_MPC,err_y_MPC,err_theta_MPC);

disp(['Quadratic Error for MPC was ',num2str(Qe_MPC)]);