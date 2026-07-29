function Z_neq_out=order_tnprod_rest(G)
N = size(G,1);
Out = G{1};
for i=1:N-2
   Out = tensor_contraction(Out,G{i+1},2*i+2,1);
end
%
Z_neq_out = tensor_contraction(Out,G{N},[1,2*N],[4,1]);
%
% Z_neq_out = reshape(Out_C,[Nway(1), prod(Nway(2:end-3)), Nway(end)]);
end