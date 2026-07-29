% mode-k unfolding of X (square bracket unfolding)
% Ik x Ik+1 ... IN I1...Ik-1
function X_sb_k=tenmat_sb(X,k)
S=size(X);
N=numel(S);
if k==1
    X_sb_k=reshape(X,[S(1), prod(S(2:end))]);
elseif k==N
     X_sb_k=reshape(X,[prod(S(1:end-1)),S(N)]);
     X_sb_k=permute( X_sb_k,[2,1]);
else
    X=reshape(X,[prod(S(1:k-1)),prod(S(k:end))]);
    X=permute(X,[2,1]);
    X_sb_k=reshape(X,[S(k),prod(S)/S(k)]);
end
%
end 