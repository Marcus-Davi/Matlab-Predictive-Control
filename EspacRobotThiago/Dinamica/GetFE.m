function [Ej,Fj] = GetFE(C,D,k)

% This function returns the values of Ej and Fj for all
% the number of pre determined interations

% Thiago Alves Lima
% Electrical Engineering Department, Federal University of Ceara

%P = [1 zeros(1,length(A))];

Ej = zeros(k,k);
Fj = zeros(k,length (D)-1);

for l = 1:k;
    v = C(1);
    Ej(l:k,l) = C(1);
    for m = 1: length (D);
        C(m) = C(m) - (v*D(m));
    end;
    C = [C(2:(length(C))) 0];
    Fj(l,:) = C(1:length (D)-1);
end;


end

