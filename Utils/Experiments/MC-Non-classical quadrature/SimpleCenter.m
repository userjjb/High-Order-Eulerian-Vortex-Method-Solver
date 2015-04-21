clear all
%Setup the test problem
del=0.01;
N=6;
[Qx,Qw] = GLquad(N);
nd=(Qx+1)/2;

[Qx2,Qw2] = GLquad(2*N);
nd2=(Qx2+1)/2;

a=0.2;
b=0.2;
dx=0.6; 
dy=0.6;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

[xx,yy]=meshgrid(nd,nd);
W=w(xx,yy);
[xx2,yy2]=meshgrid(nd2,nd2);

nn=elim(nd',nd',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
interp_w=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(W'*Lag(reshape(y,1,[]),1:N)))',size(x));

W2=interp_w(xx2,yy2);

for nx=1:N
Tx=nd(nx);
fprintf('%i ',nx)
    for ny=1:1
    Ty=0;

k=@(x,y) (y-Ty)./((x-Tx).^2+(y-Ty).^2+del^2).^(3/2);
K=k(xx,yy);
K2=k(xx2,yy2);

E(ny,nx)=integral2(@(x,y) w(x,y).*k(x,y),0,1,0,1);
I(ny,nx)=integral2(@(x,y) interp_w(x,y).*k(x,y),0,1,0,1);
Q(ny,nx)=Qw*(W.*K)*Qw'/4;
Q2(ny,nx)=Qw2*(W2.*K2)*Qw2'/4;

for Sn=1:N 
    ry=@(y) y./((y-Ty).^2+(nd(Sn)-Tx).^2+del^2).^(3/2);
    [QxJ, QwJ] = GenOrthog(2*N,ry,0,1,1);
    QxY(:,Sn)=QxJ;
    QwY(Sn,:)=QwJ;
end

WY= w(repmat(nd',2*N,1),QxY);
KY= k(repmat(nd',2*N,1),QxY);
QY(ny,nx)= sum(QwY'.*(WY.*KY))*(Qw'/2);
    end
end

difEQ=E-Q;
difEQ2=E-Q2;
difEY=E-QY;