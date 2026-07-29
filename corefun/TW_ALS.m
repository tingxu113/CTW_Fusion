function [node,Core] = TW_ALS(X,Lrank,Rrank,K)
        maxit=K;
        c=X;
        n = size(c);
        n = n(:);
        d = 2;
        node=cell(3,1);
        r=Rrank(:);
        l=Lrank(:);
        rng('default')
        node{1}=randn(r(1),n(1),l(1),r(2));
        node{2}=randn(r(2),n(2),l(2),r(3));
        node{3}=randn(r(3),n(3),l(3),r(1));
        
        Core=randn(l(1),l(2),l(3));
        Ndim=3;
        Nway=size(X);

%         TempE=(1./prod(l)).*ones(1,prod(l));
%         Core = reshape(TempE,size(Core));
        for it=1:maxit
          %% compute b
%          Girest = tnreshapemat(order_tnprod_rest(node), Ndim);
%          TempE = reshape(X,[1,prod(Nway)])/Girest;
%          Core = reshape(TempE,size(Core));
%          s=norm(Core(:));
%          Core=Core./s;
          for num = 1:Ndim
          GCrest = tenremat(circ_tnprod_rest(node,Core,num), Ndim);  % Q is the right part of the relation equation
          TempA  = tenmat_sb(X,num)/GCrest;
          node{num} = my_Fold(TempA,size(node{num}),2);
          s=norm(node{num}(:));
          node{num}=node{num}./s;
          end
         Girest = tnreshapemat(order_tnprod_rest(node), Ndim);
         TempE = reshape(X,[1,prod(Nway)])/Girest;
         Core = reshape(TempE,size(Core));
         s=norm(Core(:));
         Core=Core./s;
          
        end
%         node{3}=node{3}.*s;
        Core=Core.*s;
        t.node=node;
        t.d=d;
        t.n=n;
        t.r=r;
        return;
    end

