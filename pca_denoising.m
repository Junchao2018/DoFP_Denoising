%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% Our code is based on the following paper,So:
% Please cite the paper when you use the code:
%
% Lei Zhang, Weisheng Dong, David Zhang, Guangming Shi 
% Two-stage Image Denoising by Principal Component Analysis with Local
% Pixel Grouping, Pattern Recognition, vol. 43, issue 4, pp. 1531-1549,
% Apr. 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dI=pca_denoising(nI,v)
s=2;
[h, w]=size(nI);
%%%initial denoising
v2     =   v^2;
S      =   20;       %training block (2k+1)*(2k+1)
t      =   3;        %variable block (2t+1)*(2t+1)
nblk   =   250;
b      =   2*t+1;
b2     =   b*b;

k     =  0;
N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
r     =  [1:s:N];
r     =  [r r(end)+1:N];
c     =  [1:s:M];
c     =  [c c(end)+1:M];
X     =  zeros(b*b,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  nI(i:h-b+i,j:w-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';            
    end
end

XT       =   X';

I        =   (1:L);
I        =   reshape(I, N, M);
N1       =   length(r);
M1       =   length(c);
L        =   N1*M1;
Y        =   zeros( b2, L );

for  i  =  1 : N1
    for  j  =  1 : M1
        
        row      =   r(i);
        col      =   c(j);
        off      =   (col-1)*N + row;
        off1     =   (j-1)*N1 + i;        
        
        indc              =   LPG_new( XT, row, col, off, nblk, S, I );    
        [coe, P, V, mX]   =   getpca( X(:, indc) );
        py                =   mean(coe.^2, 2);
        px                =   max(0, py-v2);
        wei               =   px./py;
        Y(:,off1)         =   P'*(coe(:,1).*wei) + mX(:,1);
    
    end
end

% Output the processed image
dI       =  zeros(h,w);
im_wei   =  zeros(h,w);
k        =  0;
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        dI(r-1+i,c-1+j)      =  dI(r-1+i,c-1+j) + reshape( Y(k,:)', [N1 M1]);
        im_wei(r-1+i,c-1+j)  =  im_wei(r-1+i,c-1+j) + 1;
    end
end
dI        =   dI./(im_wei+eps);
end
function  indc  =  LPG_new(X, row, col, off, nv, S, I)
[N M]   =   size(I);
f2      =   size(X,2);

rmin    =   max( row-S, 1 );
rmax    =   min( row+S, N );
cmin    =   max( col-S, 1 );
cmax    =   min( col+S, M );
         
idx     =   I(rmin:rmax, cmin:cmax);
idx     =   idx(:);
B       =   X(idx, :);        
v       =   X(off, :);
        
        
dis     =   (B(:,1) - v(1)).^2;
for k = 2:f2
    dis   =  dis + (B(:,k) - v(k)).^2;
end
dis   =  dis./f2;
[val,ind]   =  sort(dis);        
indc        =  idx( ind(1:nv) );
end
function [Y, P, V, mx]=getpca(X)

%X: MxN matrix (M dimensions, N trials)
%Y: Y=P*X
%P: the transform matrix
%V: the variance vector

[M,N]=size(X);

mx=mean(X,2);
mx=repmat(mx,1,N);

X=X-mx;

CovX=X*X'/(N-1);
[P,V]=eig(CovX);

V=diag(V);
[t,ind]=sort(-V);
V=V(ind);
P=P(:,ind);
P=P';
Y=P*X;
end

