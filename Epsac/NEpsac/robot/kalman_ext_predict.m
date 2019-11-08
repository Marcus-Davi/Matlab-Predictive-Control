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
function [x_est,Pk] = kalman_ext_predict(Ts,f,Q,x,u,Pk)

%Predict
x_est = f(x,u,Ts);


x_est = reshape(x_est,length(x_est),1);


Jf = [1 0 -sin(x_est(3))*Ts; %Jacobiana do processo f(x,u)
      0 1 cos(x_est(3))*Ts;
      0 0 1];

Pk = Jf*Pk*Jf'+Q;

end


