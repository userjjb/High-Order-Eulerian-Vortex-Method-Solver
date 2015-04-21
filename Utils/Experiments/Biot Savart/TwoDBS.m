clear all
N=5;
[Qx,Qw] = GLquad(N);
nd=(Qx+1)/2;

[Qx2,Qw2] = GLquad(2*N);
nd2=(Qx2+1)/2;

a=0.3;
b=0.3;
dx=0.5; 
dy=0.5;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

[xx,yy]=meshgrid(nd,nd);
W=w(xx,yy);
[xx2,yy2]=meshgrid(nd2,nd2);

nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
interp_w=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(W'*Lag(reshape(y,1,[]),1:N)))',size(x));

W2=interp_w(xx2,yy2);

del=.5;

for nx=1:N
Tx=nd(nx);
%fprintf('%i ',nx)
    for ny=1:N
        Ty=nd(ny);

        k=@(x,y) (y-Ty)./((x-Tx).^2+(y-Ty).^2+del^2).^(3/2);
        K=k(xx,yy);
        K2=k(xx2,yy2);
        
        nyi= N-(ny-1);
        E(nyi,nx)=integral2(@(x,y) w(x,y).*k(x,y),0,1,0,1);
        I(nyi,nx)=integral2(@(x,y) interp_w(x,y).*k(x,y),0,1,0,1);
        Q(nyi,nx)=Qw*(W.*K)*Qw'/4;
        Q2(nyi,nx)=Qw2*(W2.*K2)*Qw2'/4;
    end
end

difEQ=E-Q;
difEQ2=E-Q2;