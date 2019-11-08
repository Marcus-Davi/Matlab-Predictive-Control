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
function [x_est,Kk,Pk,err] = kalman_ext_update(Ts,h,R,x,u,Pk,y)
%Update

y_measure = reshape(y,length(y),1);

Jh = eye(3); %Jacobianas da medição h(x)

yest =  h(x,u,Ts);
err = y_measure - yest;

Kk = Pk*Jh'/(Jh*Pk*Jh'+R);

x_est = x + Kk*(err);
Pk = (eye(length(Pk)) - Kk*Jh)*Pk;

end


