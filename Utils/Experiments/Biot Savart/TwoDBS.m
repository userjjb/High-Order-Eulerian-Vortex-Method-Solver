clear all
load('Iwm66BS.mat')
N=6;
f=7/36;%Size of element
del=0.0001*f;
f=f/2;

[Qx,Qw] = GLquad(N);
nd=f*(Qx+1)/1;

a=2*(f*0.4*3)^2;
b=2*(f*0.3*3)^2;
dx=f*0.6;
dy=f*0.8;
w=@(x,y) 15*exp(-((x-dx).^2/a+(y-dy).^2/b));%Center is at (dx,dy)

[xx,yy]=meshgrid(nd,nd);
W=w(xx,yy);
nn=elim(nd(1:N)',nd(1:N)',[1 3 2]);
Lag= @(x,nv) prod(bsxfun(@rdivide,bsxfun(@minus,x,nn(nv,:,:)),bsxfun(@minus,nd(nv),nn(nv,:,:))),3);
interp_w=@(x,y) reshape(diag(Lag(reshape(x,1,[]),1:N)'*(W'*Lag(reshape(y,1,[]),1:N)))',size(x));

M=N;
[Qx2,Qw2]= LGLquad(M);
nd2=f*(Qx2+1)/1;

for nx=1:M
Tx=0*f+nd2(nx); %nd3
    for ny=1:N
        Ty=0*f+nd(ny); %nd
        Rsq=@(x,y) (x-Tx).^2+(y-Ty).^2;
        %k=@(x,y) (y-Ty)./(2*pi*((x-Tx).^2+(y-Ty).^2)+delp); %RM
        %k=@(x,y) (1/(2*pi))*((y-Ty).*(Rsq(x,y)+ (2)*del^2) )./(Rsq(x,y)+del^2).^(2); %WL
        %k=@(x,y) ((y-Ty)./(2*pi*Rsq(x,y))).*(1-exp(-Rsq(x,y)/del^2)); %Gaussian
        k=@(x,y) ((y-Ty)./(2*pi*Rsq(x,y))).*(1-(1-(Rsq(x,y)/del^2)).*exp(-Rsq(x,y)/del^2)); %SG
        %k=@(x,y) ((y-Ty)./(2*pi*Rsq(x,y))).*(1-(1-(2*Rsq(x,y)/del^2+0.5*Rsq(x,y).^2/del^4)).*exp(-Rsq(x,y)/del^2)); %SG6
        %k=@(x,y) ((y-Ty)./(2*pi*Rsq(x,y))).*(1-besselj(0,sqrt(Rsq(x,y))/del)); %PS
        %k=@(x,y) ((y-Ty)./(2*pi*Rsq(x,y))).*(1- (8./(45*Rsq(x,y)/del^2)).*(4*besselj(2,4*sqrt(Rsq(x,y))/del)-5*besselj(2,2*sqrt(Rsq(x,y))/del)+besselj(2,sqrt(Rsq(x,y))/del))); %PS2
        K=k(xx,yy);
        
        nyi= N-(ny-1); %N
        I(nyi,nx)=integral2(@(x,y) interp_w(x,y).*k(x,y),0,2*f,0,2*f);%,'RelTol',1e-6,'AbsTol',1e-10);
        Q(nyi,nx)=Qw*(W.*K)*Qw'/f^-2;
        Im(nyi,nx)= sum(sum(W.*Iwm(:,:,ny,nx)))*(2*f);%*6?
    end
end

difIQ=I-Q;
%[QwS, ~,~]= StiffnessQuadModWeights(nd,nd3);

%Ltwo=sqrt(Qw*difEQ.^2*Qw3'/4)
Ltwo=sqrt(Qw*difIQ.^2*Qw2'/f^-2)*(.5/f)