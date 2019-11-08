% Argumentos 
%Ts -> tempo de amostragem
%f -> modelo estado
%h -> modelo saida
%Jf -> Jacobiana F
%Jh -> Jacobiana H
%Q -> Covarianca estado
%R -> Covarianca Saida
%x -> estado anterior
%u -> entrada 
%P -> covarianca erro
%y -> leitura
function [x_est,Kk,Pk,err] = kalman_ext(Ts,f,h,Q,R,x,u,Pk,y)

%Predict
x_est = f(x,u,Ts);


x_est = reshape(x_est,length(x_est),1);
y_measure = reshape(y,length(y),1);

Jf = [1 0 -sin(x(3))*Ts; %Jacobiana do processo f(x,u)
      0 1 cos(x(3))*Ts;
      0 0 1];

Pk = Jf*Pk*Jf'+Q;
%Update

Jh = eye(3); %Jacobianas da medição h(x)

yest =  h(x_est,u,Ts);
err = y_measure - yest;

Kk = Pk*Jh'/(Jh*Pk*Jh'+R);

x_est = x_est + Kk*(err);
Pk = (eye(length(Pk)) - Kk*Jh)*Pk;






end


