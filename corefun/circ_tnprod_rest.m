function Z_neq_out = circ_tnprod_rest(G,C,n)
G = circshift(G,-n); % arrange Z{n} to the last core, so we only need to multiply the first N-1 core
N = size(G,1);
Out = G{1};
for i=1:N-2
   Out = tensor_contraction(Out,G{i+1}, 2*i+2, 1);
end
%
m = circshift(1:N,-n);
%
Z_neq_out = tensor_contraction(Out, C, [3:2:2*N-1], m(1:N-1));
end