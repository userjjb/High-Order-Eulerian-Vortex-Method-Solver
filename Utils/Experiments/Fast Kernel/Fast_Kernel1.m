clear all
N=3;
Np=1+N;
delX=0.2;
Ex=-1:delX:1;

[Qx,Qw]=LGLquad(Np);
nd=(Qx+1)/(2/delX);

temp=bsxfun(@plus,nd,Ex(1:end-1));
xv=[reshape(temp(1:end-1,:),[],1);Ex(end)];

[xM,yM]=meshgrid(xv,xv);

del=0.1;
Ky=@(x,y) (xM-x)./((xM-x).^2+(yM-y).^2+del^2).^(3/2);

Ky11=Ky(0,0);
Ky12=Ky(nd(2),0);
Ky21=Ky(0,nd(2));
Ky13=Ky(nd(3),0);
Ky22=Ky(nd(2),nd(2));