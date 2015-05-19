clear all
N=3;
f=7.875/72;%Size of element
del=1.55*(7.875/72);

f=f/2;

[Qx,Qw] = GLquad(N);
nd=f*(Qx+1)/1;

[Qx2,Qw2] = GLquad(2*N);
nd2=f*(Qx2+1)/1;

a=2*(f*0.3)^2;
b=2*(f*0.3)^2;
dx=f*0.8; 
dy=f*0.8;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

[xx,yy]=meshgrid(nd,nd);
W=w(xx,yy);
[xx2,yy2]=meshgrid(nd2,nd2);

nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
interp_w=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(W'*Lag(reshape(y,1,[]),1:N)))',size(x));

W2=interp_w(xx2,yy2);

M=N-1;
[Qx3,Qw3]= LGLquad(M);
nd3=f*(Qx3+1)/1;


xv=0:(f*2)/10:f*2;
for nx=1:M
Tx=nd3(nx); %nd3
%fprintf('%i ',nx)
    for ny=1:N
        Ty=nd(ny); %nd
        k=@(x,y) (y-Ty)./((x-Tx).^2+(y-Ty).^2+del^2).^(3/2);
        K=k(xx,yy);
        K2=k(xx2,yy2);
        
        nyi= N-(ny-1); %N
        E(nyi,nx)=integral2(@(x,y) w(x,y).*k(x,y),0,2*f,0,2*f);
        I(nyi,nx)=integral2(@(x,y) interp_w(x,y).*k(x,y),0,2*f,0,2*f);
        Q(nyi,nx)=Qw*(W.*K)*Qw'/f^-2;
        Q2(nyi,nx)=Qw2*(W2.*K2)*Qw2'/f^-2;
    end
end

difEQ=E-Q;
difIQ=I-Q;
difEQ2=E-Q2;
%[QwS, ~,~]= StiffnessQuadModWeights(nd,nd3);

%Ltwo=sqrt(Qw*difEQ.^2*Qw3'/4)
Ltwo=sqrt(Qw*difIQ.^2*Qw3'/f^-2)*(.5/f)