function y = ajeita_G(M)
n = length(M);
y = zeros(n,n);
M2 = fliplr(M); %troca colunas
for i=1:n
    row = M2(i,:);
    row = [row zeros(1,n-i)];
    row = row(n-i+1:end);
    y(:,i) = row;
end

y = y';

end