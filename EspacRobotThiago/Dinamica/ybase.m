function [yb] = ybase(ub,ho,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[~,N] = size(ub);
[mn,nn] = size(n);

if nn == 1
    for l = 1:N-1
        n = [n n];
    end
end

for k=1:N;
    yb(:,k) = real(modelo(ub(:,k),ho,0.1*k)+n(:,k));
    ho = yb(:,k);
end

