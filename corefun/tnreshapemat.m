function Out = tnreshapemat(Grest, N)
Nway = size(Grest);
m = zeros(N,1); n = zeros(N,1);
for k=1:N
    m(k)=2*k;       n(k)=2*k-1;
end
tempG = permute(Grest,[m,n]);
Out = reshape(tempG,prod(Nway(m)),prod(Nway(n)));
%
end



