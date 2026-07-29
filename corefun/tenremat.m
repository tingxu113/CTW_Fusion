function Out = tenremat(Grest, N)
Nway = size(Grest);
%
m = [N+1,N+2,1];   n = 2:N;
%
tempG = permute(Grest,[m,n]);
Out = reshape(tempG,prod(Nway(m)),prod(Nway(n)));
end



