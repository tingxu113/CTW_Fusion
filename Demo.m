%% ================================================================
% This is the demo code for 
% Hyperspectral and Multispectral Image Fusion via Coupled Tensor Wheel Decomposition
% Ting Xu, Jin-Liang Xiao, Ya-Ru Fan, Liu Liu, and Liang-Jian Deng*
% IEEE Transactions on Geoscience and Remote Sensing, 2026.

% Contact email: tingxu@cdut.edu.cn
% If you use this code, please cite the following paper:

% @ARTICLE{11538234,
%   author={Xu, Ting and Xiao, Jin-Liang and Fan, Ya-Ru and Liu, Liu and Deng, Liang-Jian},
%   journal={IEEE Transactions on Geoscience and Remote Sensing}, 
%   title={Hyperspectral and Multispectral Image Fusion via Coupled Tensor Wheel Decomposition}, 
%   year={2026},
%   volume={64},
%   number={},
%   pages={5515114-5515114},
%   keywords={Modeling;Tensors;Heart rate;Algorithms;Superresolution;Hyperspectral imaging;Joining processes;Ranking (statistics);Aluminum;Matrices;Image fusion;proximal alternating minimization;remote sensing;tensor wheel (TW) decomposition;variational model},
%   doi={10.1109/TGRS.2026.3697926}}


% =========================================================================
%% parameters
% Please adjust the following parameters for better results
% Please adjust k at [20,80]
% Please adjust r1 at [2,11]
% Please adjust r2 at [125,400]
% Please adjust r3 at [2,11]
% Please adjust r4 at [2,4]
% =========================================================================


clear
clc
addpath(genpath('corefun'))
addpath(genpath('quality'));

load('cave_flowers.mat')
S=cave_flowers;
S=S/max(S(:));


 [M,N,L] = size(S);
 sf = 4;


% % generate HSI
 BW = ones(sf,1)/sf;
 BW1 = psf2otf(BW,[M 1]);
 BH = ones(sf,1)/sf;
 BH1 = psf2otf(BH,[N 1]);
 s0 = sf/2;
 S_w = ifft(fft(S).*repmat(BW1,1,N,L)); %blur with the width  mode
 aa = fft(permute(S_w,[2 1 3]));
 S_h = (aa.*repmat(BH1,1,M,L));
 S_h = permute(ifft(S_h),[2 1 3]);  %blur with the height mode

 HSI = S_h(s0:sf:end,s0:sf:end,:);

% % generate MSI
 F=create_F;
 S_bar = hyperConvert2D(S);
 Y=F*S_bar;

 MSI = hyperConvert3D(Y, M, N);

%%============================CTW==================================================
r1=2; r2=400; r3=5; r4=3; k=50;

par.Rrank=[r1 r2 r1];  
par.Lrank=[r3 r3 r4]; 
par.Nways=[M,N,L];
par.k=k;

t0=clock
[CTW_HSR] = TW_SR(HSI,MSI,F,BW1,BH1,sf,par,s0);
t_CTW=etime(clock,t0);

[psnr_CTW, rmse_CTW, ergas_CTW, sam_CTW, uiqi_CTW, ssim_CTW, DD_CTW, CC_CTW] = quality_assessment(CTW_HSR*255, S*255, 0, 1.0/sf);


