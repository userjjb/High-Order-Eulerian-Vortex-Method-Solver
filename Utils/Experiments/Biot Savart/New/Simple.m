clear all
%Setup the test problem
N=7;
[Qx,Qw] = GLquad(N);
nd=(Qx+1)/2;

a=0.3;
b=0.3;
dx=0.2; 
dy=0.5;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

[xx,yy]=meshgrid(nd,nd);
W=w(xx,yy);

nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
interp_w=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(W'*Lag(reshape(y,1,[]),1:N)))',size(x));

for nx=1:N
Tx=nd(nx);
fprintf('%i ',nx)
    for ny=1:3
%Target coordinates
Ty=-nd(ny);

k=@(x,y) (x-Tx)./((x-Tx).^2+(y-Ty).^2).^(3/2);
K=k(xx,yy);

E(ny,nx)=integral2(@(x,y) w(x,y).*k(x,y),0,1,0,1);
I(ny,nx)=integral2(@(x,y) interp_w(x,y).*k(x,y),0,1,0,1);
Q(ny,nx)=Qw*(W.*K)*Qw'/4;

%Calculate a spatially varying modified ortho basis
mid=0.01;
r=@(x) x./((x-Ty).^2+(mid).^2).^(3/2);
[QxSeed,QwSeed] = GenOrthog(N,r,0,1);
for Sn=1:N
    %Tensor product of spatially varying Q_r and a particular Q_r
    r=@(x) 1./((x-Tx).^2+(QxSeed(Sn)-Ty).^2).^(3/2);
    [QxJ, QwJ] = GenOrthog(N,r,0,1);
    QxM(:,Sn)=QxJ;
    QwM(Sn,:)=QwJ;
    
    %Tensor product of spatially varying Q_r and GLquad
    r=@(x) 1./((x-Tx).^2+(nd(Sn)-Ty).^2).^(3/2);
    [QxJ, QwJ] = GenOrthog(N,r,0,1);
    QxL(:,Sn)=QxJ;
    QwL(Sn,:)=QwJ;
end

QxM=QxM';
QxL=QxL';

WM= w(QxM,repmat(QxSeed',N,1)');
KM= k(QxM,repmat(QxSeed',N,1)');
QM(ny,nx)= sum(QwM'.*(WM.*KM)')*QwSeed';

WL= w(QxL,repmat(nd',N,1)');
KL= k(QxL,repmat(nd',N,1)');
QL(ny,nx)= sum(QwL'.*(WL.*KL)')*Qw'/2;
    end
end

difEQ=E-Q;
difEM=E-QM;
difEL=E-QL;