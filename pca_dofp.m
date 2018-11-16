function B=pca_dofp(I,D,par)
w = par.s;
num = par.num;
r = par.r;
[n,m]=size(I);
N=(n-w)/2+1;
M=(m-w)/2+1;
L=N*M;
X=zeros(w*w,L);
k=0;
DN=zeros(1,w*w);
for i=1:w
   for j=1:w
      k=k+1;
      X(k,:)=reshape(I(i:2:n-w+i,j:2:m-w+j),1,L);
      DN(k)=D(i,j);
   end
end
q=(L+1)/2;
Xc=X(:,q);
XC=repmat(Xc,1,L);
E=abs(X-XC);
mE=mean(E);
[val,ind]=sort(mE);
%%
X=X(:,ind(1:num));
%%
[Y, P, V, mX] =Get_PCA_Matrix(X,DN);
Y1=0*Y;
for i=1:w*w-r
   y=Y(i,:);
   p=P(i,:);
   p=p.^2;
   nv=sum(p.*DN);
   py=mean(y.^2)+0.01;
   t=max(0,py-nv);
   c=t/py;
   Y1(i,:)=c*Y(i,:);
end
%%
B=(P'*Y1+mX);
%%
B=B(:,1);
B=reshape(B,w,w);
B=B';
end

function [Y, P, V, mx]=Get_PCA_Matrix(X,D)

%X: MxN matrix (M dimensions, N trials)
%Y: Y=P*X
%P: the transform matrix
%V: the variance vector

[M,N]=size(X);
mx=mean(X,2);
mx=repmat(mx,1,N);
X=X-mx;
CovX=X*X'/(N-1);
D=diag(D);
CovX=CovX-D;
ind = find(CovX<0); 
CovX(ind) =0.0001; 
[P,V]=eig(CovX);
V=diag(V);
[t,ind]=sort(-V);
V=V(ind);
P=P(:,ind);
P=P';
Y=P*X;
end
