% ========================================================================
% PCA based denoising for DoFP image, Version 1.0
% Copyright(c) 2017  Junchao Zhang, Haibo Luo, Rongguang Liang, 
% Wei Zhou, Bin Hui and Zheng Chang
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%
%----------------------------------------------------------------------
% This is an demo of "PCA-based denoising method for division of focal plane
% polarimeters"
% 
% Please cite the following paper when you use it:
%
% Junchao Zhang, Haibo Luo, Rongguang Liang, Wei Zhou, Bin Hui and 
% Zheng Chang, "PCA-based denoising method for division of focal plane 
% polarimeters," Optics Express 25, 2391-2400 (2017).
% ========================================================================

clc;
close all;
clear all;
%% Load Image
img1=double(imread('./Test image/0.bmp'));
img2=double(imread('./Test image/45.bmp'));
img3=double(imread('./Test image/90.bmp'));
img4=double(imread('./Test image/135.bmp'));
[rows,cols]=size(img1);
%% Calculate ground-truth of Stokes parameters
S0_ori=(img1+img2+img3+img4)*0.5;
S1_ori=img1-img3;
S2_ori=img2-img4;
DOLP_ori=sqrt(S1_ori.^2+S2_ori.^2)./S0_ori;
I_ori=zeros(rows,cols);
I_ori(1:2:end,1:2:end)=img1(1:2:end,1:2:end);
I_ori(1:2:end,2:2:end)=img2(1:2:end,2:2:end);
I_ori(2:2:end,1:2:end)=img4(2:2:end,1:2:end);
I_ori(2:2:end,2:2:end)=img3(2:2:end,2:2:end);
%% Interpolation
[B0, B45, B90, B135]=interpolation(I_ori);
S0_interp=(B0+B90+B45+B135)*0.5;
S1_interp=B0-B90;
S2_interp=B45-B135;
DOLP_interp=sqrt(S1_interp.^2+S2_interp.^2)./S0_interp;
%% Crop image
img1=img1(233:end,89:436);
img2=img2(233:end,89:436);
img3=img3(233:end,89:436);
img4=img4(233:end,89:436);
I_ori=I_ori(233:end,89:436);
S0_ori=S0_ori(233:end,89:436);
DOLP_ori=DOLP_ori(233:end,89:436);
S0_interp=S0_interp(233:end,89:436);
DOLP_interp=DOLP_interp(233:end,89:436);
%% Add noise
v1=1;
v2=2;
v3=2;
v4=1;
[img1,u1]=addnoise(img1,v1,17);
[img2,u2]=addnoise(img2,v2,18);
[img3,u3]=addnoise(img3,v4,19);
[img4,u4]=addnoise(img4,v3,20);
%% noise DOFP data
[rn,cn]=size(img1);
In=zeros(rn,cn);
In(1:2:end,1:2:end)=img1(1:2:end,1:2:end);
In(1:2:end,2:2:end)=img2(1:2:end,2:2:end);
In(2:2:end,1:2:end)=img4(2:2:end,1:2:end);
In(2:2:end,2:2:end)=img3(2:2:end,2:2:end);
%% PCA-based ÂË³ýËæ»úÔëÉù
I=In;
[n,m]=size(I);
%% parameters
s=6;% block size
k=34;% training block size,(k-s)/2 should be an even integer
k2=k/2;
c=1.1;
par.s = s;
par.num = 50;% the number of the first similarity blcoks. 
par.r = 8; % the rest rows set to zero.
%%%type1
D_type1=zeros(s,s);
D_type1(1:2:s,1:2:s)=c*v1;
D_type1(1:2:s,2:2:s)=c*v2;
D_type1(2:2:s,1:2:s)=c*v3;
D_type1(2:2:s,2:2:s)=c*v4;
D_type1=D_type1.^2;
%%%type2
D_type2=zeros(s,s);
D_type2(1:2:s,1:2:s)=c*v2;
D_type2(1:2:s,2:2:s)=c*v1;
D_type2(2:2:s,1:2:s)=c*v4;
D_type2(2:2:s,2:2:s)=c*v3;
D_type2=D_type2.^2;
%%%type3
D_type3=zeros(s,s);
D_type3(1:2:s,1:2:s)=c*v3;
D_type3(1:2:s,2:2:s)=c*v4;
D_type3(2:2:s,1:2:s)=c*v1;
D_type3(2:2:s,2:2:s)=c*v2;
D_type3=D_type3.^2;
%%%type4
D_type4=zeros(s,s);
D_type4(1:2:s,1:2:s)=c*v4;
D_type4(1:2:s,2:2:s)=c*v3;
D_type4(2:2:s,1:2:s)=c*v2;
D_type4(2:2:s,2:2:s)=c*v1;
D_type4=D_type4.^2;
im_wei=zeros(n,m);
dItmp=zeros(n,m);
%% denoising
tic;
for i=1:1:n-k
   for j=1:1:m-k
      Block=I(i:i+k-1,j:j+k-1);%k by k block
      if(mod(i,2)~=0&&mod(j,2)~=0)
          D=D_type1;
      else if(mod(i,2)~=0&&mod(j,2)==0)
              D=D_type2;
          else if(mod(i,2)==0&&mod(j,2)~=0)
                  D=D_type3;
              else
                  D=D_type4;
              end
          end
      end
      dB=pca_dofp(Block,D,par);%pca denoising
      dItmp(i-1+k2:i+k2,j-1+k2:j+k2)=dItmp(i-1+k2:i+k2,j-1+k2:j+k2)+dB(3:4,3:4);
      im_wei(i-1+k2:i+k2,j-1+k2:j+k2)=im_wei(i-1+k2:i+k2,j-1+k2:j+k2)+1;
   end
end
dI=dItmp./(im_wei+eps);
toc;
%% Calculate PSNR
% I
e=I_ori-In;
e=e(30+1:end-30,30+1:end-30);
me=mean(mean(e.^2));
sI=10*log10(255^2/me);
% S0
[B0, B45, B90, B135]=interpolation(In);
S0_n=(B0+B90+B45+B135)*0.5;
S1_n=B0-B90;
S2_n=B45-B135;
DOLP_n=sqrt(S1_n.^2+S2_n.^2)./S0_n;
e=S0_ori-S0_n;%e=S0_ori-S0_n;
e=e(30+1:end-30,30+1:end-30);
me=mean(mean(e.^2));
sS0=10*log10(510^2/me);
% DOLP
e=DOLP_ori-DOLP_n;%e=DOLP_ori-B_DOLP;
e=e(30+1:end-30,30+1:end-30);
me=mean(mean(e.^2));
sDolp=10*log10(1^2/me);
disp(['Befor denoising<I  S0  DoLP>--is:',num2str(sI),'   ',num2str(sS0),'   ',num2str(sDolp)]);
%%
%
e=I_ori-dI;
e=e(30+1:end-30,30+1:end-30);
me=mean(mean(e.^2));
sI=10*log10(255^2/me);
% S0
[B0, B45, B90, B135]=interpolation(dI);
S0_dn=(B0+B90+B45+B135)*0.5;
S1_dn=B0-B90;
S2_dn=B45-B135;
DOLP_dn=sqrt(S1_dn.^2+S2_dn.^2)./S0_dn;
e=S0_ori-S0_dn;%e=S0_ori-S0_dn;
e=e(30+1:end-30,30+1:end-30);
me=mean(mean(e.^2));
sS0=10*log10(510^2/me);
%DOLP
e=DOLP_ori-DOLP_dn;%e=DOLP_ori-DOLP_dn;
e=e(30+1:end-30,30+1:end-30);
me=mean(mean(e.^2));
sDolp=10*log10(1^2/me);
disp(['After denoising<I  S0  DoLP>--is:',num2str(sI),'   ',num2str(sS0),'   ',num2str(sDolp)]);
%% Show results
figure;imshow(S0_ori(30+1:end-30,30+1:end-30),[55,220]);
figure;imshow(DOLP_ori(30+1:end-30,30+1:end-30),[0.01,0.26]);

figure;imshow(S0_interp(30+1:end-30,30+1:end-30),[55,220]);
figure;imshow(DOLP_interp(30+1:end-30,30+1:end-30),[0.01,0.26]);

figure;imshow(S0_n(30+1:end-30,30+1:end-30),[55,220]);
figure;imshow(DOLP_n(30+1:end-30,30+1:end-30),[0.01,0.26]);

figure;imshow(S0_dn(30+1:end-30,30+1:end-30),[55,220]);
figure;imshow(DOLP_dn(30+1:end-30,30+1:end-30),[0.01,0.26]);
%% Residual noise remove
K=30;
dif=dI(K+1:end-K,K+1:end-K)-In(K+1:end-K,K+1:end-K);
dif1=dif(1:2:end,1:2:end);
dif2=dif(1:2:end,2:2:end);
dif3=dif(2:2:end,1:2:end);
dif4=dif(2:2:end,2:2:end);

v1=v1^2-(mean(mean(dif1.^2)));
v2=v2^2-(mean(mean(dif2.^2)));
v3=v3^2-(mean(mean(dif3.^2)));
v4=v4^2-(mean(mean(dif4.^2)));

v1=sqrt(abs(v1));
v2=sqrt(abs(v2));
v3=sqrt(abs(v3));
v4=sqrt(abs(v4));

sigma_S0=sqrt((v1^2+v2^2+v3^2+v4^2)*0.25);
sigma_S1=sqrt(v1^2+v3^2);
sigma_S2=sqrt(v2^2+v4^2);
%%
c=1.1;
S0_sn=pca_denoising(S0_dn,sigma_S0*c);
S1_sn=pca_denoising(S1_dn,sigma_S1*c);
S2_sn=pca_denoising(S2_dn,sigma_S2*c);
%%
DOLP_sn=sqrt(S1_sn.^2+S2_sn.^2)./S0_sn;
figure;imshow(S0_sn(30+1:end-30,30+1:end-30),[55,220]);
figure;imshow(DOLP_sn(30+1:end-30,30+1:end-30),[0.01,0.26]);
%
e=S0_ori-S0_sn;%e=S0_ori-S0_dn;
e=e(50+1:end-50,50+1:end-50);
me=mean(mean(e.^2));
sS0=10*log10(510^2/me);
%DOLP
e=DOLP_ori-DOLP_sn;%e=DOLP_ori-DOLP_dn;
e=e(50+1:end-50,50+1:end-50);
me=mean(mean(e.^2));
sDolp=10*log10(1^2/me);
disp(['After residual noise removing<~  S0  DoLP>--is:',num2str(sS0),'   ',num2str(sDolp)]);