clear all
%Setup the test problem
N=6;
[Qx,Qw] = GLquad(N);
nd=(Qx+1)/2;

a=0.2;
b=0.2;
dx=0.5; 
dy=0.5;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

[xx,yy]=meshgrid(nd,nd);
W=w(xx,yy);

nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
interp_w=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(W'*Lag(reshape(y,1,[]),1:N)))',size(x));

for nx=1:N
Tx=nd(nx)
    for ny=1:N
%Target coordinates
Ty=-nd(ny);

k=@(x,y) (x-Tx)./((x-Tx).^2+(y-Ty).^2).^(3/2);
K=k(xx,yy);

E(ny,nx)=integral2(@(x,y) w(x,y).*k(x,y),0,1,0,1);
I(ny,nx)=integral2(@(x,y) interp_w(x,y).*k(x,y),0,1,0,1);
Q(ny,nx)=Qw*(W.*K)*Qw'/4;

%Calculate a spatially varying modified ortho basis
mid=0;
r=@(x) 1./((x-Tx).^2+(mid-Ty).^2).^(3/2);
[QxSeed,QwSeed] = GenOrthog(N,r,0,1);
for Sn=1:N
    %Tensor product of spatially varying Q_r and a particular Q_r
    r=@(x) 1./((x-Ty).^2+(QxSeed(Sn)-Tx).^2).^(3/2);
    [QxJ, QwJ] = GenOrthog(N,r,0,1);
    QxM(:,Sn)=QxJ;
    QwM(Sn,:)=QwJ;
    
    %Tensor product of spatially varying Q_r and GLquad
    r=@(x) 1./((x-Ty).^2+(nd(Sn)-Tx).^2).^(3/2);
    [QxJ, QwJ] = GenOrthog(N,r,0,1);
    QxL(:,Sn)=QxJ;
    QwL(Sn,:)=QwJ;
end
mid2=0;
r=@(x) 1./((x-Tx).^2+(mid2-Ty).^2).^(3/2);
[QxT,QwT] = GenOrthog(N,r,0,1);

WM= w(repmat(QxSeed',N,1),QxM);
KM= k(repmat(QxSeed',N,1),QxM);
QM(ny,nx)= sum(QwM'.*(WM.*KM))*QwSeed';

WL= w(repmat(nd',N,1),QxL);
KL= k(repmat(nd',N,1),QxL);
QL(ny,nx)= sum(QwL'.*(WL.*KL))*Qw'/2;

WT= w(repmat(QxT',N,1),repmat(QxT,1,N));
KT= k(repmat(QxT',N,1),repmat(QxT,1,N));
QT(ny,nx)= QwT*(WT.*KT)*QwT';

    end
end

difEQ=E-Q;
difET=E-QT;
difEM=E-QM;
difEL=E-QL;