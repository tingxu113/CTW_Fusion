 function [HR_HSI,rmse1] = TW_SR(HSI,MSI,T,BW,BH,downsampling_scale,par,s0)
mu = 1;
Lrank = par.Lrank;
Rrank = par.Rrank;
R=[Rrank;Lrank];                  
Nways=par.Nways;
%%  simulate LR-HSI
Y_h_bar = hyperConvert2D(HSI);
HSI1 = tenmat_sb(HSI,1);
HSI2 = tenmat_sb(HSI,2);
HSI3 = tenmat_sb(HSI,3);
%%  simulate HR-MSI
MSI1 = tenmat_sb(MSI,1);
MSI2 = tenmat_sb(MSI,2);
MSI3 = tenmat_sb(MSI,3);

%% inilization D1 D2 D3
Ndim=3;
G = cell(Ndim,1);
[G,~] = TW_ALS(MSI,Lrank,Rrank,10);

%%
Core=rand([Lrank(1),Lrank(2),Lrank(3)]);
Y_h_bar = hyperConvert2D(HSI);
N3 = My_vca(Y_h_bar,Lrank(3)*Rrank(3)*Rrank(1));
N3=reshape(N3,[Nways(3) Lrank(3) Rrank(3) Rrank(1)]);
G{3}=permute(N3,[3 1 2 4]);

%%

D1 = G{1}; D2 = G{2}; %D3 = G{3};
D_1 = ifft(fft(my_Unfold(D1,size(D1),2)).*repmat(BW,[1 Rrank(1)*Rrank(2)*Lrank(1)]));

D_1 = my_Fold(D_1(s0:downsampling_scale:end,:),[Rrank(1),size(HSI,1),Lrank(1),Rrank(2)],2);

D_2 = ifft(fft(my_Unfold(D2,size(D2),2)).*repmat(BH,[1 Rrank(2)*Rrank(3)*Lrank(2)]));
D_2 = my_Fold(D_2(s0:downsampling_scale:end,:),[Rrank(2),size(HSI,2),Lrank(2),Rrank(3)],2);


D3=my_Unfold(G{3},size(G{3}),2);
D_3 = my_Fold(T*D3,[Rrank(3),size(MSI,3),Lrank(3),Rrank(1)],2);

D3=G{3};
D11 = cell(3,1);
D22 = cell(3,1);
D11{1} = D_1;
D11{2} = D_2;
D11{3} = D3;
D22{1} = D1;
D22{2} = D2;
D22{3} = D_3;
%% iteration
rmse1=zeros(80,1);
kk=par.k;

for i=1:kk

%%  update D1
QH = tenremat(circ_tnprod_rest(D11,Core,1), 3);
QM = tenremat(circ_tnprod_rest(D22,Core,1), 3);
D1_M =my_Unfold(D1,size(G{1}),2);
[ D1 ] = CG1_TW( MSI1, D1_M, QM,HSI1,downsampling_scale,s0,BW,QH,mu);
D_1 = real(ifft(fft(D1).*repmat(BW,[1 Rrank(1)*Rrank(2)*Lrank(1)])));
D1 =my_Fold(D1,size(G{1}),2);
D_1 =my_Fold(D_1(s0:downsampling_scale:end,:),[Rrank(1),size(HSI,1),Lrank(1),Rrank(2)],2);%张量形式的P1G1
D11{1} = D_1;D22{1} = D1;

%%  update D2
QH = tenremat(circ_tnprod_rest(D11,Core,2), 3);
QM = tenremat(circ_tnprod_rest(D22,Core,2), 3);
D2_M =my_Unfold(D2,size(G{2}),2);
[ D2 ] = CG1_TW( MSI2, D2_M, QM,HSI2,downsampling_scale,s0,BH,QH,mu );
D_2 = real(ifft(fft(D2).*repmat(BH,[1 Rrank(2)*Rrank(3)*Lrank(2)])));
D2 =my_Fold(D2,size(G{2}),2); 
D_2 =my_Fold(D_2(s0:downsampling_scale:end,:),[Rrank(2),size(HSI,2),Lrank(2),Rrank(3)],2);%张量形式的P1G1
D11{2} = D_2;D22{2} = D2;


%%  update D3

QH = tenremat(circ_tnprod_rest(D11,Core,3), 3);
QM = tenremat(circ_tnprod_rest(D22,Core,3), 3);
D3_M =my_Unfold(D3,size(D3),2);
[ D3 ] = CG2_TW( HSI3, D3_M, QH,MSI3,T,QM,mu );


D_3 =my_Fold(T*D3,[Rrank(3),size(MSI,3),Lrank(3),Rrank(1)],2);%张量形式的P1G1
D3 =  my_Fold(D3,[Rrank(3),size(HSI,3),Lrank(3),Rrank(1)],2);
D11{3} = D3;D22{3} = D_3;


%%  update Core 
QH = tnreshapemat(order_tnprod_rest(D11), 3);
QM = tnreshapemat(order_tnprod_rest(D22), 3);
TempC = reshape(HSI,[1,prod(size(HSI))])*QH'+reshape(MSI,[1,prod(size(MSI))])*QM'+mu*reshape(Core,[1,numel(Core)]);

TempD = (QH*QH')+(QM*QM')+(mu)*eye(size(QH,1),size(QH,1));
TempE = TempC*pinv(TempD);
Core = reshape(TempE,size(Core));

end
 GG = cell(3,1); GG{1} = D1; GG{2} = D2; GG{3} = D3;
 Girest=tnreshapemat(order_tnprod_rest(GG), 3);
 TempE=reshape(Core,[1,numel(Core)]);
 HR_HSI= reshape(TempE*Girest,Nways);
end

