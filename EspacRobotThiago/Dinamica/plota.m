legends{1,1} = strcat('Reference robot');
legends{2,1} = strcat('EPSAC-Sub');


figure(1)
plot(xr(1:end-5*N),yr(1:end-5*N),'k','LineWidth',2); hold on;
plot(x_EPSAC,y_EPSAC,'b:','LineWidth',1.5); hold on;
legend(legends,'Location','northwest','Interpreter','latex');
title('Trajectory');
ylabel('$Y$ (m)','Interpreter','latex');
xlabel('$X$ (m)','Interpreter','latex');
axis tight;
grid on;


figure(2)
plot(t,tr(1:end-5*N),'k','LineWidth',2); hold on;
plot(t,theta_EPSAC,'b:','LineWidth',1.5); hold on;
legend(legends,'Location','southeast','Interpreter','latex');
title('Orientation');
ylabel('$\theta$ (rad)','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
axis tight;
grid on;
%

legends{1,1} = strcat('EPSAC-Sub');

figure(3)
subplot(3,1,1);
plot(t,err_x_EPSAC,'b:','LineWidth',1.5); hold on;
legend('EPSAC-Sub','Kanklar','Interpreter','latex');
ylabel('$(x_r-x)$ (m)','Interpreter','latex');
title('Pose erros');
axis tight;
grid on;
subplot(3,1,2);
plot(t,err_y_EPSAC,'b:','LineWidth',1.5); hold on;
ylabel('$(y_r-y)$ (m)','Interpreter','latex');
axis tight;
grid on;
subplot(3,1,3);
plot(t,err_theta_EPSAC,'b:','LineWidth',1.5); hold on;
axis tight;
grid on;
ylabel('$(\theta_r-\theta$ (rad))','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
%
legends{1,1} = strcat('Reference robot');
legends{2,1} = strcat('EPSAC-Sub');

figure(4)
subplot(2,1,1);

plot(t,vr(1:end-5*N),'k','LineWidth',2); hold on;
v_EPSAC = [0;v_EPSAC(1:end-1)];
stairs(t,v_EPSAC,'b:','LineWidth',1.5); hold on;
legend(legends,'Location','southeast','Interpreter','latex');
title('Control signals','Interpreter','latex');
ylabel('$\nu$ (m/s)','Interpreter','latex');
axis tight;
grid on;

subplot(2,1,2);
plot(t,wr(1:end-5*N),'k','LineWidth',2); hold on;
w_EPSAC = [0;w_EPSAC(1:end-1)];
stairs(t,w_EPSAC,'b:','LineWidth',1.5); hold on;
xlabel('Time (s)','Interpreter','latex');
ylabel('$\omega$ (rad/s)','Interpreter','latex');
axis tight;
grid on;

